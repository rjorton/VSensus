package vsensus;

import java.io.*;
import java.util.*;
import java.text.DecimalFormat;

public class VSENSUS {

    public static DecimalFormat df4 = new DecimalFormat("#.####");
    public static int totalAmbCount=0;
    
    public static void main(String[] args) {

        System.out.println("VSENSUS viral consensus caller started...\n");
        
        //tidy up - nucleotide or base in messages
        
        String mpileupFileArg="", seqHeaderArg="";
        int minCovArg=5;
        int qFilArg=0;
        int qOffsetArg=33;
        double ambiFreqArg=0.4;
        double indelFreqArg=0.49;
        boolean ambiguitiesArg=true, insertionsArg=true, refPadArg=false, gapPadArg=false, noIndelsArg=false;        
        
        String usage="Required:\njava -jar VSENSUS.jar mpileup.txt";
        usage+="\n\nOptional:\njava -jar VSENSUS.jar mpileup.txt -c=minCov -q=baseQual -a=ambiguityFreq -i=indelFreq -v -o -m -r -g -n";
        
        usage+="\n\nRequired arguments:";
        usage+="\nmpileup.txt : mpileup file created via samtools from a BAM file";
        usage+="\n\nOptional arguments:";
        usage+="\n-c=minCov : [integer] the minimum Coverage at a genome position with which to call a consensus nucleotide: default >= "+minCovArg;
        usage+="\n-q=baseQual : [integer] the minimum nucleotide Quality for inclusion: default >= "+qFilArg;
        usage+="\n-a=ambiguityFreq : [double] the secondary nucleotide frequency needed to trigger Ambiguity codes in the consensus: default >= "+ambiFreqArg;
        usage+="\n-i=indelFreq : [double] frequency at which to call Indels into the consensus: default = "+indelFreqArg;
        usage+="\n-v : aVoid using ambiguities if possible [not possible if two or more nucleotides are tied in count and quality]";
        usage+="\n-o : use a phred 64 Offest rather than 33";
        usage+="\n-m : use the Mode insertion of the Mode insertion length, rather than the plain mode insertion";
        usage+="\n-r : at sites below the minimum cov use the Reference base rather than an N in the consensus";
        usage+="\n-g : at sites below the minimum cov use Gap (-) rather than an N in the consensus [overrides -p]";
        usage+="\n-n : do Not call indels at all [not recommended - but sometimes useful]";
        usage+="\n-h=seqHeader : set fasta consensus sequence Header \n";
 
        if(args.length>0) {
                    
            mpileupFileArg=args[0];
            
            for(int i=1;i<args.length;i++) {
                if(args[i].indexOf("-c=")==0) {
                    //minCov >=
                    try {
                        minCovArg=Integer.parseInt(args[i].substring(3, args[i].length()));
                    } catch (NumberFormatException e) {
                        System.out.println("\nError - unrecognised integer for -c: "+args[i]);
                        System.out.println(usage);
                        System.exit(1);
                    }
                    
                    if(minCovArg<1) {
                        System.out.println("Warning -c=minCov set below 1 - setting to 1");
                        minCovArg=1;
                    }
                }
                else if(args[i].indexOf("-q=")==0) {
                    //min baseQual >=
                    try {
                        qFilArg=Integer.parseInt(args[i].substring(3, args[i].length()));
                    } catch (NumberFormatException e) {
                        System.out.println("\nError - unrecognised integer for -q: "+args[i]);
                        System.out.println(usage);
                        System.exit(1);
                    }
                    
                    if(qFilArg<0) {
                        System.out.println("Warning -q=baseQual set below 0 - setting to 0");
                        qFilArg=0;
                    }
                }
                else if(args[i].indexOf("-a=")==0) {
                    //frequency at which to consider ambiguity >=
                    try {
                        ambiFreqArg=Double.parseDouble(args[i].substring(3, args[i].length()));
                    } catch (NumberFormatException e) {
                        System.out.println("\nError - unrecognised double for -a: "+args[i]);
                        System.out.println(usage);
                        System.exit(1);
                    }
                    
                    if(ambiFreqArg>1 | ambiFreqArg<0) {
                        System.out.println("\nError -a=ambiguityFreq set to greater than 1 or less than 0: "+args[i]);
                        System.out.println(usage);
                        System.exit(1);
                    }
                    
                    if(ambiFreqArg>0.5) {
                        System.out.println("Warning -a=ambiguityFreq set greater than 0.5 1 - setting to 0.5");
                        ambiFreqArg=0.5;
                    }
                }
                else if(args[i].indexOf("-i=")==0) {
                    //indelFreq
                    try {
                        indelFreqArg=Double.parseDouble(args[i].substring(3, args[i].length()));
                    } catch (NumberFormatException e) {
                        System.out.println("Error - unrecognised double for -a: "+args[i]);
                        System.out.println(usage);
                        System.exit(1);
                    }
                    
                    if(indelFreqArg>1 | indelFreqArg<0) {
                        System.out.println("\nError -i=indelFreq set to greater than 1 or less than 0: "+args[i]);
                        System.out.println(usage);
                        System.exit(1);
                    }
                }
                else if(args[i].indexOf("-h=")==0) {
                    seqHeaderArg=args[i].substring(3, args[i].length());
                }
                else if(args[i].indexOf("-o")==0) {
                    //qOffset set to 64, default 33
                    qOffsetArg=64;
                }
                else if(args[i].indexOf("-v")==0) {
                    //avoid ambiguites - so set ambiguities to false
                    ambiguitiesArg=false;
                }
                else if(args[i].indexOf("-m")==0) {
                    //insertion mode - use mode of mode length
                    insertionsArg=false;
                }
                else if(args[i].indexOf("-r")==0) {
                    //pad consensus with the refBase at positions below minCov
                    refPadArg=true;
                }
                else if(args[i].indexOf("-g")==0) {
                    //pad consensus with - at positions below minCov
                    gapPadArg=true;
                    refPadArg=false;
                }
                else if(args[i].indexOf("-n")==0) {
                    //do not call indels
                    noIndelsArg=true;
                }
                else {
                    System.out.println("\nError - unrecognised argument: "+args[i]);
                    System.out.println("Outputting help file:");
                    System.out.println("\n"+usage);
                    System.exit(1);
                }
            }
            
            createConsensus(mpileupFileArg, minCovArg, qFilArg, ambiFreqArg, qOffsetArg, ambiguitiesArg, insertionsArg, refPadArg, indelFreqArg, gapPadArg, noIndelsArg, seqHeaderArg);
        }
        else {
           System.out.println("\nError - no arguments provided");
           System.out.println("Outputting help file:");
           System.out.println("\n"+usage);
           System.exit(1);
        }
        
        System.out.println("\n...VSENSUS consensus caller finished");
    }
    

    public static void createConsensus(String mpileupFile, int minCov, int qFil, double ambiFreq, int qOffset, boolean ambiguities, boolean modeInsertions, boolean refPad, double indelFreq, boolean gapPad, boolean noIndels, String seqHeader) {
        
        File inFile = new File(mpileupFile);
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));
            input.close();
        }    
        catch (IOException ex) {
            System.out.println("Error - mpileup file does not exist: "+mpileupFile);
            System.exit(1);
        }
                
        System.out.println("VSENSUS settings...");
        System.out.println("-Input BAM mpileup file: "+mpileupFile);
        System.out.println("-Minimum coverage >= "+minCov);
        System.out.println("-Minimum base quality >= "+qFil);
        System.out.println("-Phred quality score offset = "+qOffset);
        
        if(ambiguities) {
            System.out.println("-Amiguity code base frequency >= "+ambiFreq);
        }
        else {
            System.out.println("-Avoid ambiguities in the consensus sequence [but unavoidable if same base count and quality]");
        }
        
        if(noIndels) {
            System.out.println("-Indel calling has been disabled");
        }
        else {
            System.out.println("-Indel frequency threshold >= "+indelFreq);
            
            if(modeInsertions) {
                System.out.println("-Insertion selection will utilise the overall mode insertion");
            }
            else {
                System.out.println("-Insertion selection will first utilise the mode insertion of the mode insertion length");
            }
        }
        
        if(gapPad) {
            System.out.println("-Using a - [rather than an N or reference base] at zero coverage positions");
        }
        else if(refPad) {
            System.out.println("-Using the reference base [rather than an N or -] at zero coverage positions");
        }
        else {
            System.out.println("-Using an N [rather than the reference base or -] at zero coverage positions");
        }
        
        String diversityFile, consensusFile;
        
        if(mpileupFile.lastIndexOf(".")>0) {
            diversityFile=mpileupFile.substring(0, mpileupFile.lastIndexOf("."))+"_vsensus.txt";
            consensusFile=mpileupFile.substring(0, mpileupFile.lastIndexOf("."))+"_vsensus.fasta";
        }
        else {
            diversityFile=mpileupFile+"_vsensus.txt";
            consensusFile=mpileupFile+"_vsensus.fasta";
        }
        
        System.out.println("-Output VSENSUS diversity file = "+diversityFile);
        System.out.println("-Output VSENSUS consensus file = "+consensusFile);
        
        if(seqHeader.equals("")) {
            String stub=mpileupFile;
            int loc=mpileupFile.indexOf("_mpileup.txt");
            if(loc>0)
                stub=mpileupFile.substring(0, loc);
            seqHeader="VSENSUS-"+stub;
        }
        System.out.println("-Output fasta sequence header = >"+seqHeader+"_referenceSequenceName");
        
        System.out.println("...VSENSUS settings\n");

        System.out.println("Calling consenus");
        
        String conSeq="";
        double indelWarnFreq=0.25;
        
        try {
            FileWriter fstream = new FileWriter(diversityFile);
            BufferedWriter out = new BufferedWriter(fstream);

            String header="Chr"
                    +"\tPos"
                    + "\tRef"
                    + "\tPileCov"
                    + "\tBaseCov"
                    + "\tQFilt"
                    + "\tFwdCov"
                    + "\tRevCov"
                    + "\tRefN"
                    + "\tNonRefN"
                    + "\tMismatch"
                    + "\tMismatchNonZero"
                    + "\tTopMutation"
                    + "\tEntropy"
                    + "\tAvQ"
                    + "\tAvQ-P"
                    + "\tA"
                    + "\tC"
                    + "\tG"
                    + "\tT"
                    + "\tN"
                    + "\tDelProp"
                    + "\tDelStart"
                    + "\tDels"
                    + "\tMeanDelLen"
                    + "\tModeDelLen"
                    + "\tModeDelFreq"
                    + "\tInsProp"
                    + "\tInsStart"
                    + "\tMeanInsLen"
                    + "\tModeInsLen"
                    + "\tModeInsLenFreq"
                    + "\tModeIns"
                    + "\tModeInsFreq"
                    //+ "\tAltModeIns"
                    //+ "\tAltModeInsFreq"
                    + "\tReadsStart"
                    + "\tReadsEnd"
                    + "\tAfwd"
                    + "\tCfwd"
                    + "\tGfwd"
                    + "\tTfwd"
                    + "\tNfwd"
                    + "\tArev"
                    + "\tCrev"
                    + "\tGrev"
                    + "\tTrev"
                    + "\tNrev"
                    + "\tAqual"
                    + "\tCqual"
                    + "\tGqual"
                    + "\tTqual"
                    + "\tNqual";
            
            out.write(header+"\n"); 

            inFile = new File(mpileupFile);
            try {
                BufferedReader input =  new BufferedReader(new FileReader(inFile));

                char nts[]=new char[10];
                nts[0]='A';
                nts[1]='C';
                nts[2]='G';
                nts[3]='T';
                nts[4]='N';
                nts[5]='a';
                nts[6]='c';
                nts[7]='g';
                nts[8]='t';
                nts[9]='n';

                int skip=0, runningPos=0, chrCount=0;
                
                String prevChrName="";
                boolean prevChange=false;
                
                double entropyFreq=0, entropyTotal=0;
                int entropyCount=0, entropyZeroCount=0;
                
                int totalMutCount=0, totalDelCount=0, totalInsCount=0, totalAlertCount=0, totalIndelWarnCount=0;
                int aboveMinCov=0, belowMinCov=0, totalPositions=0;
                
                try {
                    String line = null;

                    while ((line = input.readLine()) != null) {
                        
                        int pos=0, pileCov=0, refBasePos=0, removedBases=0;
                        int fwdCov=0, revCov=0;
                        int nonRefCount=0, refCount=0;
                        int mapCount=0, mapSum=0, readEnds=0;
                        int baseCounts[]=new int[10];
                        int finalCounts[]=new int[5];
                        double finalFreqs[]=new double[5];
                        
                        int qualCount=0, qualSum=0;
                        double qualAv=0, qualPSum=0, qualPAv=0;
                        String qualString="";
                        int baseQuals[]=new int[10];
                        double baseQualPs[]=new double[10];
                        double baseQualsAv[]=new double[10];
                        double finalQualsAv[]=new double[5];
                        double finalQualPsAv[]=new double[5];
                        
                        int insCount=0, insSum=0;
                        double insMeanLen=0, insProp=0;
                        int delCount=0, delSum=0, delStar=0;
                        double delMeanLen=0, delProp=0;
                        
                        char refBase='?';
                        
                        String chrName="", outLine="";
                        boolean change=false;
                        
                        HashMap<Integer, Integer> delLenHash = new HashMap<>();
                        HashMap<Integer, Integer> insLenHash = new HashMap<>();
                        HashMap<String, Integer> insHash = new HashMap<>();

                        String[] splits=line.split("\t");

                        for(int i=0;i<splits.length;i++) {

                            if(i==0) {
                                chrName=splits[i];
                            }
                            else if(i==1) {
                                pos=Integer.parseInt(splits[i]);
                                totalPositions++;
                                
                                if(pos==1 | !prevChrName.equals(chrName)) {
                                    System.out.println("Chromosome = "+chrName);
                                    
                                    if(chrCount>0)
                                        conSeq+="\n";
                                    
                                    conSeq+=">"+seqHeader+"_"+chrName+"\n";
                                    chrCount++;
                                }                         
                            }
                            else if(i==2) {
                                refBase=splits[i].toUpperCase().charAt(0);
                                refBasePos=-1;

                                for(int j=0;j<5;j++) {
                                    if(refBase==nts[j]) {
                                        refBasePos=j;
                                        break;
                                    }
                                }

                                if(refBasePos<0) {
                                    System.out.println("Error - ambiguous character in reference at pos "+pos+" = "+refBase+" - setting to N");
                                    refBasePos=4;
                                }
                            }

                            else if(i==3) {
                                pileCov=Integer.parseInt(splits[i]);
                            }

                            else if(i==4) {
                                if(pileCov==0 & splits[4].equals("*")) {
                                    splits[4]="";
                                    splits[5]="";
                                }
                                //splits[4] is the sequence - but load the qualString [5] as well
                                qualString=splits[5];
                                int qualPos=-1;
                                char thisBase='?';

                                for(int j=0;j<splits[i].length();j++) {

                                    thisBase='?';

                                    boolean found=false;

                                    //find the ACGT acgt mutations
                                    for(int k=0;k<nts.length;k++) {
                                        if(splits[i].charAt(j)==nts[k]) {
                                            qualPos++;

                                            if((int)qualString.charAt(qualPos)-qOffset>=qFil) {
                                                baseCounts[k]++;
                                                thisBase=nts[k];
                                                found=true;
                                            }
                                            else {
                                                removedBases++;
                                            }
                                            
                                            break;
                                        }
                                    }

                                    if(found) {
                                        //mutation ACGTN acgtn - already found above
                                    }
                                    else if(splits[i].charAt(j)=='.') {
                                        qualPos++;
                                        if((int)qualString.charAt(qualPos)-qOffset>=qFil) {
                                            baseCounts[refBasePos]++;
                                            thisBase=nts[refBasePos];
                                            found=true;
                                        }
                                        else {
                                            removedBases++;
                                        }
                                    }
                                    else if(splits[i].charAt(j)==',') {
                                        qualPos++;
                                        if((int)qualString.charAt(qualPos)-qOffset>=qFil) {
                                            baseCounts[refBasePos+5]++;
                                            thisBase=nts[refBasePos+5];
                                            found=true;
                                        }
                                        else {
                                            removedBases++;
                                        }
                                    }
                                    else if(splits[i].charAt(j)=='^') {
                                        //mapping quality ^Phred
                                        //have to increment the seqfield j along another one to account for Q score itself
                                        mapCount++;
                                        mapSum+=(int)splits[i].charAt(j+1)-qOffset;
                                        j++;
                                    }
                                    else if(splits[i].charAt(j)=='+') {
                                        insCount++;
                                        int insStringLen=1;
                                        boolean test=true;
                                        
                                        while(test) {
                                            char c=splits[i].charAt(j+insStringLen);

                                            if(c>='0' & c<='9') {
                                                insStringLen++;
                                            }
                                            else {
                                                int insLen=Integer.parseInt(splits[i].substring(j+1, j+insStringLen));
                                                insSum+=insLen;
                                                String insString=splits[i].substring(j+insStringLen, j+insStringLen+insLen).toUpperCase();
                                                j+=insLen+insStringLen-1;
                                                
                                                if(!insHash.containsKey(insString)) {
                                                    insHash.put(insString,1);
                                                }
                                                else {
                                                    insHash.put(insString, insHash.get(insString)+1);
                                                }
                                                
                                                if(!insLenHash.containsKey(insLen)) {
                                                    insLenHash.put(insLen,1);
                                                }
                                                else {
                                                    insLenHash.put(insLen, insLenHash.get(insLen)+1);
                                                }
                                                
                                                test=false;
                                            }
                                        }
                                    }
                                    else if(splits[i].charAt(j)=='-') {
                                        delCount++;
                                        int delStringLen=1;
                                        
                                        boolean test=true;

                                        while(test) {
                                            char c=splits[i].charAt(j+delStringLen);

                                            if(c>='0' & c<='9') {
                                                delStringLen++;
                                            }
                                            else {
                                                int delLen=Integer.parseInt(splits[i].substring(j+1, j+delStringLen));
                                                delSum+=delLen;
                                                j+=delLen+delStringLen-1;
                                                
                                                if(!delLenHash.containsKey(delLen)) {
                                                    delLenHash.put(delLen,1);
                                                }
                                                else {
                                                    delLenHash.put(delLen, delLenHash.get(delLen)+1);
                                                }
                                                
                                                test=false;
                                            }
                                        }
                                    }
                                    else if(splits[i].charAt(j)=='*') {
                                        //mpileup gives a placeholder quality score of deletion - so increment qualPos aling one
                                        delStar++;
                                        qualPos++;
                                        thisBase='.';
                                        found=true;
                                    }
                                    else if(splits[i].charAt(j)=='$') {
                                        readEnds++;
                                    }
                                    else {
                                        //turning off for Q filter
                                        //i.e. a q filtered mutated base will not be found nor will it be any of the above
                                        //System.out.println("ERROR - char not recognised "+splits[i].charAt(j));
                                    }
                                
                                    if(found) {//if it is real reference base (not an read start, end, insertion etc)
                                        if(thisBase!='.') {//ignore the delStars - for some reason they have qualities here
                                            int thisQual=(int)qualString.charAt(qualPos)-qOffset;
                                            qualCount++;
                                            qualSum+=thisQual;
                                            double qProb=Math.pow((double)10,(double)-thisQual/(double)10);
                                            qualPSum+=qProb;
                                            
                                            for(int k=0;k<nts.length;k++) {
                                                if(thisBase==nts[k]) {
                                                    baseQuals[k]+=thisQual;
                                                    baseQualPs[k]+=qProb;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    
                                }//end of seq/base length (splits[4]) loop
                            }//end of i==4 mpileup field
                        }//end of splits loop to read in data

                        if(insCount>0)
                            insMeanLen=(double)insSum/(double)insCount;
                        if(delCount>0)
                            delMeanLen=(double)delSum/(double)delCount;
                        if(qualCount>0)
                            qualAv=(double)qualSum/(double)qualCount;
                        if(qualPSum>0)
                             qualPAv=qualPSum/(double)qualCount;

                        int baseCov=0;
                        for(int i=0;i<finalCounts.length;i++) {

                            finalCounts[i]=baseCounts[i]+baseCounts[i+5];
                            baseCov+=finalCounts[i];

                            if(i==refBasePos) {
                                refCount=finalCounts[i];
                            }
                            else {
                                nonRefCount+=finalCounts[i];
                            }
                            
                            fwdCov+=baseCounts[i];
                            revCov+=baseCounts[i+5];
                        }
                        
                        for(int i=0;i<finalCounts.length;i++) {
                            finalFreqs[i]=(double)finalCounts[i]/(double)baseCov;
                        }
                        
                        for(int i=0;i<baseCounts.length;i++) {
                            if(baseCounts[i]>0) {
                                baseQualsAv[i]=(double)baseQuals[i]/(double)baseCounts[i];
                            }
                        }
                        
                        for(int i=0;i<finalQualsAv.length;i++) {
                            if(baseCounts[i]>0 | baseCounts[i+5]>0) {
                                finalQualsAv[i]=(double)(baseQuals[i]+baseQuals[i+5])/(double)(baseCounts[i]+baseCounts[i+5]);
                                finalQualPsAv[i]=(double)(baseQualPs[i]+baseQualPs[i+5])/(double)(baseCounts[i]+baseCounts[i+5]);
                            }
                        }

                        double entropy=0;
                        for(int i=0;i<finalFreqs.length-1;i++) {
                            if(finalFreqs[i]>0 & finalFreqs[i]<1) {
                                entropy-=finalFreqs[i]*Math.log(finalFreqs[i]);
                            }
                        }
                        
                        //only include in the overall entropy totals if meets minCov threshold
                        if(baseCov>=minCov) {
                            entropyTotal+=entropy;
                            entropyCount++;

                            if(entropy==0) {
                                entropyZeroCount++;
                            }   
                        }

                        int maxBaseCount=0;
                        for(int i=0;i<finalCounts.length;i++) {
                            if(finalCounts[i]>maxBaseCount) {
                                maxBaseCount=finalCounts[i];
                            }
                        }
                        //its max as using phred scores not the Ps (which would be minimum)
                        double maxBaseQual=0;
                        for(int i=0;i<finalCounts.length;i++) {
                            if(finalCounts[i]==maxBaseCount) {
                                if(finalQualsAv[i]>maxBaseQual) {
                                    maxBaseQual=finalQualsAv[i];
                                }
                            }
                        }
                        
                        //Above we found the maximumBaseCount
                        //And the maximum base (Phred) quality of the bases that are max
                        //Now we get all the bases that have the maximumBaseCount - as can be more than one
                        //If ambiguities is false - we will try and split it based on quality - but will not force it beyond that
                        //i.e. even with avoid-ambiguities there may be some
                        String ambiBases="";
                        int ambiCount=0;
                        if(maxBaseCount>0 & baseCov>=minCov) {
                            for(int i=0;i<finalCounts.length;i++) {
                                if(finalCounts[i]==maxBaseCount) {
                                    if(!ambiguities) {
                                        if(finalQualsAv[i]==maxBaseQual) {
                                            ambiCount++;
                                            ambiBases+=nts[i];
                                        }
                                    }
                                    else {
                                        ambiCount++;
                                        ambiBases+=nts[i];
                                    }
                                }
                            }
                        }
                        char conBase=getAmbiguityCode(ambiBases);
                        
                        if(!ambiguities & ambiCount>1) {
                            System.out.println("Pos "+pos+" is ambiguous [could not be resolved based on quality] = "+conBase+" [A:"+df4.format(finalFreqs[0])+" C:"+df4.format(finalFreqs[1])+" G:"+df4.format(finalFreqs[2])+" T:"+df4.format(finalFreqs[3])+" N:"+df4.format(finalFreqs[4])+"] cov="+baseCov);
                        }
                                    
                        if(ambiguities & baseCov>=minCov) {
                            ambiCount=0;
                            ambiBases="";
                            
                            for(int i=0;i<finalFreqs.length;i++) {
                                if(finalFreqs[i]>=ambiFreq) {
                                    ambiCount++;
                                    ambiBases+=nts[i];
                                }
                            }

                            if(ambiCount==0) {
                                for(int i=0;i<finalFreqs.length;i++) {
                                    if(finalFreqs[i]>=0.25) {
                                        ambiCount++;
                                        ambiBases+=nts[i];
                                    }
                                }
                                ambiBases=""+getAmbiguityCode(ambiBases);
                                System.out.println("Pos "+pos+" - warning - ambiguous - no bases were observed above the ambiguity frequency threshold "+ambiFreq+ ", so ambiFreq forced to 0.25 for this position = "+ambiBases+" [A:"+df4.format(finalFreqs[0])+" C:"+df4.format(finalFreqs[1])+" G:"+df4.format(finalFreqs[2])+" T:"+df4.format(finalFreqs[3])+" N:"+df4.format(finalFreqs[4])+"] cov="+baseCov);
                                //ambiBases=""+conBase;
                            }
                            else if(ambiCount==1) {
                                //do nothing - 1 nucleotide is winner - will be set to conBase at end
                            }
                            else {
                                ambiBases=""+getAmbiguityCode(ambiBases);
                                System.out.println("Pos "+pos+" is ambiguous = "+ambiBases+" [A:"+df4.format(finalFreqs[0])+" C:"+df4.format(finalFreqs[1])+" G:"+df4.format(finalFreqs[2])+" T:"+df4.format(finalFreqs[3])+" N:"+df4.format(finalFreqs[4])+"] cov="+baseCov);
                            }
                            
                            conBase=ambiBases.charAt(0);
                        }
                        
                        double mismatch=0;
          
                        if(baseCov>0)
                            mismatch=(double)nonRefCount/(double)baseCov;
                        
                        int maxNonRefBase=0;
                        for(int i=0;i<finalCounts.length;i++) {
                            if(i==refBasePos)
                                continue;
                            if(finalCounts[i]>maxNonRefBase)
                                maxNonRefBase=finalCounts[i];
                        }
                        double maxNonRefMismatch=0;
                        if(baseCov>0)
                            maxNonRefMismatch=(double)maxNonRefBase/(double)baseCov;

                        //1 in a million default
                        double nonZeroMismatch=mismatch;
                        if(baseCov<1 | nonZeroMismatch==0)
                            nonZeroMismatch=0.000001;
                        
                        if(baseCov!=qualCount) {
                            System.out.println("Pos "+pos+" - warning - possible error at position "+pos+ " baseCov does not equal qualCount: "+baseCov+"/"+qualCount);
                        }
                        
                        int baseDelCov=0;
                        for(int i=0;i<5;i++) {
                            baseDelCov+=baseCounts[i]+baseCounts[i+5];
                        }
   
                        baseDelCov+=delStar;

                        if(baseDelCov+removedBases!=pileCov) {
                            System.out.println("Pos "+pos+" - warning - possible error at position "+pos+ " baseDelCov does not equal mpileupCov : "+(baseDelCov+removedBases)+"/"+pileCov);
                        }
      
                        //we use the total insFreq (insCount/baseCov) to determine whether to call an insertion
                        //default is to then identify the mode length of insertion, then the mode of those
                        //alternative is to find the pure mode insertion
                        //in all cases largest insertion rules in event of a tie, and first one encountered if same length
                        
                        //we get the most common insertion length - default to the large one if a tie
                        int modeInsLen=0;
                        int modeInsLenFreq=0;
                        for (Map.Entry<Integer, Integer> entry : insLenHash.entrySet()) {
                            int hashInsLen=entry.getKey();
                            int hashInsFreq=entry.getValue();
                            
                            if(hashInsFreq>modeInsLenFreq | (hashInsFreq==modeInsLenFreq & hashInsLen>modeInsLen)) {
                                modeInsLenFreq=hashInsFreq;
                                modeInsLen=hashInsLen;
                            }
                        }
                        
                        //then find the most common insertion of the most common length
                        //if a tie - just use the first one - otherwise this will go on and on
                        int modeInsFreq=0;
                        String modeIns="";
                        for (Map.Entry<String, Integer> entry : insHash.entrySet()) {
                            String hIns=entry.getKey();
                            int hashInsFreq=entry.getValue();
                            
                            if(hIns.length()==modeInsLen) {
                                if(hashInsFreq>modeInsFreq) {
                                    modeInsFreq=hashInsFreq;
                                    modeIns=hIns;
                                }
                            }
                        }
                        
                        //find the true most common insertion
                        //if tied use the largest
                        int trueModeInsFreq=0;
                        String trueModeIns="";      
                        for (Map.Entry<String, Integer> entry : insHash.entrySet()) {
                            String hashIns=entry.getKey();
                            int hashInsFreq=entry.getValue();
                            
                            if(hashInsFreq>trueModeInsFreq | (hashInsFreq==trueModeInsFreq & hashIns.length()>trueModeIns.length())) {
                                trueModeInsFreq=hashInsFreq;
                                trueModeIns=hashIns;
                            }
                        }

                        if(modeInsertions) {
                            modeIns=trueModeIns;
                            modeInsFreq=trueModeInsFreq;
                        }
                        
                        //deletions are simpler, find most common, default to largest if tied
                        int modeDelFreq=0;
                        int modeDelLen=0;
                        for (Map.Entry<Integer, Integer> entry : delLenHash.entrySet()) {
                            int hashDelLen=entry.getKey();
                            int hashDelFreq=entry.getValue();
                            
                            if(hashDelFreq>modeDelFreq | (hashDelFreq==modeDelFreq & hashDelLen>modeDelLen)) {
                                modeDelFreq=hashDelFreq;
                                modeDelLen=hashDelLen;
                            }
                        }
                        
                        insProp=(double)insCount/(double)baseCov;
                        delProp=(double)delCount/(double)baseCov;
                        
                        outLine="";
                        outLine+=chrName
                                +"\t"+pos
                                +"\t"+refBase
                                +"\t"+pileCov
                                +"\t"+baseCov
                                +"\t"+removedBases
                                +"\t"+fwdCov
                                +"\t"+revCov
                                +"\t"+refCount
                                +"\t"+nonRefCount
                                +"\t"+mismatch
                                +"\t"+nonZeroMismatch
                                +"\t"+maxNonRefMismatch
                                +"\t"+entropy
                                +"\t"+qualAv
                                +"\t"+qualPAv;
                        
                        for(int i=0;i<5;i++) {
                            outLine+="\t"+finalCounts[i];
                        }
                        
                        outLine+="\t"+delProp
                                +"\t"+delCount
                                +"\t"+delStar
                                +"\t"+delMeanLen
                                +"\t"+modeDelLen
                                +"\t"+modeDelFreq
                                +"\t"+insProp
                                +"\t"+insCount
                                +"\t"+insMeanLen
                                +"\t"+modeInsLen
                                +"\t"+modeInsLenFreq
                                +"\t"+modeIns
                                +"\t"+modeInsFreq
                                //+"\t"+trueModeIns
                                //+"\t"+trueModeInsFreq
                                +"\t"+mapCount
                                +"\t"+readEnds;
                        
                        for(int i=0;i<baseCounts.length;i++) {
                            outLine+="\t"+baseCounts[i];
                        }
                        for(int i=0;i<finalQualsAv.length;i++) {
                            outLine+="\t"+finalQualsAv[i];
                        }

                        out.write(outLine+"\n");

                        //deletion handling - if "within" a deletion don't add anything at this site - skip it
                        if(skip>0) {
                            skip--;
                            runningPos++;
                        }
                        else {
                            //old, but left in incase not used -a on mplieup
                            //catch gaps in coverage at the start and midway - pad with Ns
                            //wont handle gaps in coverage at 3' end
                            while(runningPos<(pos-1)) {
                                conSeq+="N";
                                runningPos++; 
                            }

                            if(baseCov>=minCov) {
                                aboveMinCov++;
                                
                                if(conBase==refBase) {
                                    conSeq+=refBase;
                                }
                                else {
                                   conSeq+=conBase;
                                   totalMutCount++;
                                   System.out.println("Pos "+pos+" mut "+refBase+" > "+conBase+" [A:"+df4.format(finalFreqs[0])+" C:"+df4.format(finalFreqs[1])+" G:"+df4.format(finalFreqs[2])+" T:"+df4.format(finalFreqs[3])+" N:"+df4.format(finalFreqs[4])+"] cov="+baseCov);
   
                                   if(prevChange) {
                                       totalAlertCount++;
                                       System.out.println("Alert - neighbouring sites with mutations/indels");
                                   }
                                   change=true;
                                }
                            }
                            else {
                                belowMinCov++;
                                
                                if(gapPad)
                                    conSeq+="-";
                                 else if(refPad)
                                    conSeq+=refBase;
                                else
                                    conSeq+="N";
                            }

                            if(baseCov>=minCov) {
                                if(noIndels) {
                                    if(insProp>=indelFreq) {
                                        System.out.println("Pos "+pos+" - warning - you have turned indels off - and there are insertions at this position freq = "+insProp);
                                        System.out.println("Pos "+pos+" - warning - ignoring indels could lead to incorrect consensus bases at or next to this position");
                                    }
                                    if(delProp>=indelFreq) {
                                        System.out.println("Pos "+pos+" - warning - you have turned indels off - and there are deletions at this position freq = "+delProp);
                                        System.out.println("Pos "+pos+" - warning - ignoring indels could lead to incorrect consensus bases at or next to this position");
                                    }
                                }
                                else {
                                    //cant have both an insertion and a deletion, ins first if not del
                                    if(insProp>=indelFreq) {
                                        if(trueModeIns.length()!=modeInsLen) {
                                            System.out.println("Pos "+pos+" - warning - there are more insertions of length "+modeInsLen+" than there are of the most common insertion "+trueModeIns+" length "+trueModeIns.length());
                                            System.out.println("Pos "+pos+" the insertion "+modeIns+" of the mode length "+modeInsLen+" has been selected");
                                        }
                                        conSeq+=modeIns;
                                        totalInsCount++;
                                        System.out.println("Pos "+pos+" ins "+modeInsLen+" "+modeIns+" freq="+df4.format(insProp)+" cov="+baseCov);
                                        
                                        if(prevChange) {
                                            totalAlertCount++;
                                            System.out.println("Alert - neighbouring sites with mutations/indels");
                                        }
                                        
                                        change=true;
                                    } 
                                    else if(delProp>=indelFreq) {
                                        skip=modeDelLen;
                                        totalDelCount++;
                                        System.out.println("Pos "+pos+" del "+modeDelLen+" freq="+df4.format(delProp)+" cov="+baseCov);
                                        
                                        if(prevChange) {
                                            totalAlertCount++;
                                            System.out.println("Alert - neighbouring sites with mutations/indels");
                                        }
                                        change=true;
                                    }

                                    if(delProp<indelFreq & delProp>=indelWarnFreq) {
                                        totalIndelWarnCount++;
                                        System.out.println("Pos "+pos+" - warning - there is a high proportion of [sub-consensus level] deletions: "+delProp+", mode length = "+modeDelLen+", cov="+baseCov);
                                    }

                                    if(insProp<indelFreq & insProp>=indelWarnFreq) {
                                        totalIndelWarnCount++;
                                        System.out.println("Pos "+pos+" - warning - there is a high proportion of [sub-consensus level] insertions: "+insProp+", mode insertion = "+modeIns+", cov="+baseCov);
                                    }
                                }
                            }
                            
                            //keep running pos ticking over with pos
                            runningPos++;
                        }//end of skip/else
                        
                        prevChrName=chrName;
                        prevChange=change;
                    }//end of line while
  
                    System.out.println("\nSummary");
                    System.out.println(chrCount+" = Total segments/sequences");
                    System.out.println(totalPositions+" = Total positions");
                    System.out.println(aboveMinCov+" = Total positions where coverage >= "+minCov);
                    System.out.println(belowMinCov+" = Total positions where coverage < "+minCov);
                    System.out.println(totalMutCount+" = Total mutations [includes ambiguities]");
                    System.out.println(totalAmbCount+" = Total ambiguitites");
                    System.out.println(totalInsCount+" = Total insertions");
                    System.out.println(totalDelCount+" = Total deletions");
                    System.out.println(totalIndelWarnCount+" = Total indel warnings [high freq but not consensus]");
                    System.out.println(totalAlertCount+" = Total alerts [mutations/indels at neighbouring positions]");
                    System.out.println("\nEntropy");
                    System.out.println("[Total Sites, Zero Sites, Total Entropy, Average Entropy] = ["+entropyCount+", "+entropyZeroCount+", "+entropyTotal+", "+entropyTotal/(double)entropyCount+"]");
                }
                finally {
                    input.close();
                }
            }

            catch (IOException ex) {
                ex.printStackTrace();
            }
          
            out.close();
        }

        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
         
        //Output the consneus sequence
        try {
            FileWriter fstream = new FileWriter(consensusFile);
            BufferedWriter out = new BufferedWriter(fstream);

            out.write(conSeq+"\n");
            out.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    public static char getAmbiguityCode(String ambiBases) {
        if(ambiBases.length()==0) {
            ambiBases="N";
        }
        else if(ambiBases.length()==1) {
            if(!ambiBases.equals("A") & !ambiBases.equals("C") & !ambiBases.equals("G") & !ambiBases.equals("T")) {
                ambiBases="N";
                totalAmbCount++;
            }
        }
        else if(ambiBases.length()==2) {
            if((ambiBases.charAt(0)=='A' & ambiBases.charAt(1)=='C'))
                ambiBases="M";
            else if((ambiBases.charAt(0)=='A' & ambiBases.charAt(1)=='G'))
                ambiBases="R";
            else if((ambiBases.charAt(0)=='A' & ambiBases.charAt(1)=='T'))
                ambiBases="W";
            else if((ambiBases.charAt(0)=='C' & ambiBases.charAt(1)=='G'))
                ambiBases="S";
            else if((ambiBases.charAt(0)=='C' & ambiBases.charAt(1)=='T'))
                ambiBases="Y";
            else if((ambiBases.charAt(0)=='G' & ambiBases.charAt(1)=='T'))
                ambiBases="K";
            else {
                ambiBases="N";
            }
            totalAmbCount++;
        }
        else if(ambiBases.length()==3) {
            if((ambiBases.charAt(0)=='A' & ambiBases.charAt(1)=='C' & ambiBases.charAt(2)=='G'))
                ambiBases="V";
            else if((ambiBases.charAt(0)=='A' & ambiBases.charAt(1)=='C' & ambiBases.charAt(2)=='T'))
                ambiBases="H";
            else if((ambiBases.charAt(0)=='A' & ambiBases.charAt(1)=='G' & ambiBases.charAt(2)=='T'))
                ambiBases="D";
            else if((ambiBases.charAt(0)=='C' & ambiBases.charAt(1)=='G' & ambiBases.charAt(2)=='T'))
                ambiBases="B";
            else {
                ambiBases="N";
            }
            totalAmbCount++;
        }
        else {
            ambiBases="N";
            totalAmbCount++;
        }
        
        return ambiBases.charAt(0);
    }
   
}
