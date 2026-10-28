#!/bin/bash
# Launch all jobs with dependencies
# WARNING! All scripts must be edited prior to submitting and not changed until run!
# HINT: Give meaningful job names in submission scripts to distinguish them

# Global variables
DEPENDS="--dependency=afterok:"
SCRIPTPATH="./00_scripts/slurm_code"

echo "WARNING! Make sure you have properly set up to config file in config/ folder"
echo " "

option=$1
if [ -z "$option" ] ; then
    echo "Error no option provided ! I don't know what to do"
    Help
    exit 2
else
    echo -e  "\n***will run workflow with option : $option*****\n\n"
fi

#getting parameters in config file:
source ./config/config

## checking option 
if [ "$option" = 1 ] || [ "$option" = 6 ] ;
then
   echo "test"
   #checking parameters in config files:
   if [ -z "$haplotype2" ] && [ -n "$haplotype1" ] && [ -n "$genome1" ]; then
       echo -e "only $genome1 for $haploype1 provided\n"
       #assumes only genome1 and haplotype1 are provided: 
       genome="$genome1"
       folder="haplo1"
       haplo="$haplotype1"
       if [ ! -d "$folder" ] ; then mkdir $folder ; fi
       if [ ! -d "$folder"/03_genome ] ; then mkdir $folder/03_genome ; fi
       if [ "$rnaseq" == YES ] ; then 
           p1=$(sbatch "$SCRIPTPATH"/01_submit_trimmomatic.sh                                              | awk '{print $4}')
           p2=$(sbatch "$DEPENDS"$p1 "$SCRIPTPATH"/02_submit_gmap.sh "$haplo" "$folder" "$genome"          | awk '{print $4}')
           p3=$(sbatch "$DEPENDS"$p2 "$SCRIPTPATH"/03_gsnap_array.sh "$haplo" "$folder"           	       | awk '{print $4}')
           if [ "$annotateTE" == YES ] ; then
               p4=$(sbatch "$SCRIPTPATH"/04_submit_repeatmodeler.sh -g "$genome" -s "$haplo" -f "$folder"  | awk '{print $4}')
               p5=$(sbatch "$DEPENDS"$p4:$p3 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome" -s "$haplo" -o "$folder" -r YES | awk '{print $4}')
               p5=$(sbatch "$DEPENDS"$p5 "$SCRIPTPATH"/05_braker_RNAseq.sh   -g "$genome" -s "$haplo" -o "$folder" -r YES | awk '{print $4}')
               p6=$(sbatch "$DEPENDS"$p5 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome" -s "$haplo" -o "$folder" -r NO  | awk '{print $4}')
           else
               p4=$(sbatch "$DEPENDS"$p3 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome" -s "$haplo" -o "$folder" -r YES | awk '{print $4}')
               p5=$(sbatch "$DEPENDS"$p4 "$SCRIPTPATH"/05_braker_rnaseq.sh   -g "$genome" -s "$haplo" -o "$folder" -r YES | awk '{print $4}')
               p6=$(sbatch "$DEPENDS"$p4 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome" -s "$haplo" -o "$folder" -r NO  | awk '{print $4}')
           fi
	   p7=$(sbatch "$DEPENDS"$p5:$p6 "$SCRIPTPATH"/06_reshape_braker_output.sh -s "$haplo" -g "$genome" -f "$folder" -r YES|awk '{print $4'} ) 

       elif [ "$rnaseq" == NO ]; then
           if [ "$annotateTE" == YES ] ; then
                p4=$(sbatch "$SCRIPTPATH"/04_submit_repeatmodeler.sh -g "$genome" -s "$haplo" -f "$folder"  | awk '{print $4}')
		p5=$(sbatch "$DEPENDS"$p4 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome" -s "$haplo" -o "$folder" -r NO | awk '{print $4}')
    	        p6=$(sbatch "$DEPENDS"$p5 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome" -s "$haplo" -o "$folder" -r NO  | awk '{print $4}')
            elif [ "$annotateTE" == NO ]; then
		p5=$(sbatch "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome" -s "$haplo" -o "$folder" -r YES | awk '{print $4}')
    	        p6=$(sbatch "$DEPENDS"$p5 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome" -s "$haplo" -o "$folder" -r NO  | awk '{print $4}')
       	    fi
	    p7=$(sbatch "$DEPENDS"$p6 "$SCRIPTPATH"/06_reshape_braker_output.sh -s $haplo -g $genome -f "$folder" -r NO |awk '{print $4'} ) 
       fi
  #checking parameters in config files:
  elif [ -n "$haplotype2" ] && [ -n "$haplotype1" ] && [ -n "$genome1" ] &&  [ -n "$genome2" ] ; then
    echo -e "$genome1 for $haploype1 provided\n"
    echo -e "$genome2 for $haploype2 provided\n"

	#assumes both genome1, haplotype1, genome2, haplotype2 are provided:
	folder1="haplo1"
	folder2="haplo2"
	if [ ! -d "$folder1" ] ; then mkdir $folder1 ; fi
	if [ ! -d "$folder2" ] ; then mkdir $folder2 ; fi
           p1=$(sbatch "$SCRIPTPATH"/01_submit_trimmomatic.sh                                              | awk '{print $4}')

       if [ "$rnaseq" == YES ] ; then 
	   echo -e "rnaseq provided\n"
           p2=$(sbatch "$DEPENDS"$p1 "$SCRIPTPATH"/02_submit_gmap.sh "$haplotype1" "$folder1" "$genome1"       | awk '{print $4}')
           p3=$(sbatch "$DEPENDS"$p2 "$SCRIPTPATH"/03_gsnap_array.sh "$haplotype1" "$folder1"           	   | awk '{print $4}')
	       #genome2:
           p4=$(sbatch "$DEPENDS"$p1 "$SCRIPTPATH"/02_submit_gmap.sh "$haplotype2" "$folder2" "$genome2"       | awk '{print $4}')
           p5=$(sbatch "$DEPENDS"$p4 "$SCRIPTPATH"/03_gsnap_array.sh "$haplotype2" "$folder2"           	   | awk '{print $4}')

           if [ "$annotateTE" == YES ] ; then
	           echo -e "TE annotation requested\n"	
               p6=$(sbatch "$SCRIPTPATH"/04_submit_repeatmodeler.sh -g "$genome1" -s "$haplotype1" -f "$folder1"  | awk '{print $4}') 
  	       p7=$(sbatch "$DEPENDS"$p6 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome1" -s "$haplotype1" -o "$folder1" -r YES | awk '{print $4}')
               p8=$(sbatch "$DEPENDS"$p7:$p3 "$SCRIPTPATH"/05_braker_RNAseq.sh   -g "$genome1" -s "$haplotype1" -o "$folder1" -r YES | awk '{print $4}')
               p9=$(sbatch "$DEPENDS"$p7:$p3 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome1" -s "$haplotype1" -o "$folder1" -r NO  | awk '{print $4}')
	           #genome2:
               p9=$(sbatch "$SCRIPTPATH"/04_submit_repeatmodeler.sh -g "$genome2" -s "$haplotype2" -f "$folder2"  | awk '{print $4}')
  	       p10=$(sbatch "$DEPENDS"$p9 "$SCRIPTPATH"/05_prepare_data_braker.sh  -g "$genome2" -s "$haplotype2" -o "$folder2" -r YES | awk '{print $4}')
               p11=$(sbatch "$DEPENDS"$p10:$p5 "$SCRIPTPATH"/05_braker_RNAseq.sh   -g "$genome2" -s "$haplotype2" -o "$folder2" -r YES | awk '{print $4}')
               p12=$(sbatch "$DEPENDS"$p10:$p5 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome2" -s "$haplotype2" -o "$folder2" -r NO  | awk '{print $4}')

           else
  	       p7=$(sbatch "$DEPENDS"$p6 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome1" -s "$haplotype1" -o "$folder1" -r YES | awk '{print $4}')
               p8=$(sbatch "$DEPENDS"$p7 "$SCRIPTPATH"/05_braker_RNAseq.sh   	 -g "$genome1" -s "$haplotype1" -o "$folder1" -r YES | awk '{print $4}')
               p9=$(sbatch "$DEPENDS"$p7 "$SCRIPTPATH"/05_braker_db_array.sh 	 -g "$genome1" -s "$haplotype1" -o "$folder1" -r NO  | awk '{print $4}')
	       #genome2:
  	       p10=$(sbatch "$DEPENDS"$p5 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome2" -s "$haplotype2" -o "$folder2" -r YES | awk '{print $4}')
               p11=$(sbatch "$DEPENDS"$p10 "$SCRIPTPATH"/05_braker_RNAseq.sh     -g "$genome2" -s "$haplotype2" -o "$folder2" -r YES | awk '{print $4}')
               p12=$(sbatch "$DEPENDS"$p10 "$SCRIPTPATH"/05_braker_db_array.sh    -g "$genome2" -s "$haplotype2" -o "$folder2" -r NO  | awk '{print $4}')
           fi
	   p12=$(sbatch "$DEPENDS"$p8:$p9 "$SCRIPTPATH"/06_reshape_braker_output.sh   -s $haplotype1 -g $genome1 -f "$folder1" -r YES |awk '{print $4'} ) 
	   p14=$(sbatch "$DEPENDS"$p11:$p12 "$SCRIPTPATH"/06_reshape_braker_output.sh -s $haplotype2 -g $genome2 -f "$folder2" -r YES |awk '{print $4'} ) 

       elif [ "$rnaseq" == NO ]; then
	       echo "no rnaseq provided"
           if [ "$annotateTE" == YES ] ; then
                p6=$(sbatch "$SCRIPTPATH"/04_submit_repeatmodeler.sh -g "$genome1" -s "$haplotype1" -f "$folder1"  | awk '{print $4}')
   	        p7=$(sbatch "$DEPENDS"$p6 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome1" -s "$haplotype1" -o "$folder1" -r NO | awk '{print $4}')
    	        p8=$(sbatch "$DEPENDS"$p7 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome1" -s "$haplotype1" -o "$folder1" -r NO  | awk '{print $4}')
		#genome2:
		p9=$(sbatch "$SCRIPTPATH"/04_submit_repeatmodeler.sh -g "$genome2" -s "$haplotype2" -f "$folder2"  | awk '{print $4}')
   	        p10=$(sbatch "$DEPENDS"$p9 "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome2" -s "$haplotype2" -o "$folder2" -r NO | awk '{print $4}')
    	        p11=$(sbatch "$DEPENDS"$p10 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome2" -s "$haplotype2" -o "$folder2" -r NO  | awk '{print $4}')

            elif [ "$annotateTE" == NO ] ; then
   	        p7=$(sbatch "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome1" -s "$haplotype1" -o "$folder1" -r NO | awk '{print $4}')
    	        p8=$(sbatch "$DEPENDS"$p7 "$SCRIPTPATH"/05_braker_db_array.sh  -g "$genome1" -s "$haplotype1" -o "$folder1" -r NO  | awk '{print $4}')
	        #genome2:  
		p10=$(sbatch "$SCRIPTPATH"/05_prepare_data_braker.sh -g "$genome2" -s "$haplotype2" -o "$folder2" -r NO | awk '{print $4}')
    	        p11=$(sbatch "$DEPENDS"$p10 "$SCRIPTPATH"/05_braker_db_array.sh -g "$genome2" -s "$haplotype2" -o "$folder2" -r NO  | awk '{print $4}')

	   fi
		#reshape:
	   p12=$(sbatch "$DEPENDS"$p8 "$SCRIPTPATH"/06_reshape_braker_output.sh  -s $haplotype1 -g $genome1 -f "$folder1" -r NO |awk '{print $4'} ) 
	   p14=$(sbatch "$DEPENDS"$p11 "$SCRIPTPATH"/06_reshape_braker_output.sh -s $haplotype2 -g $genome2 -f "$folder2" -r NO |awk '{print $4'} ) 
      fi
   fi
fi
if [ "$option" = 1 ] ;
then
    echo -e "option number $option is required \n\n" 
    echo "run data preparation for these option:"    
    ##NOT TESTED YET:
    if [ -z "$haplotype2" ] && [ -n "$haplotype1" ] && [ -n "$genome1" ]; then
       echo -e "only $genome1 for $haploype1 provided\n"
       p15=$(sbatch "$DEPENDS"$p7 "$SCRIPTPATH"/07_dataprep_opt3_opt4_opt5.sh |awk '{print $4}' )
    elif [ -n "$haplotype2" ] && [ -n "$haplotype1" ] && [ -n "$genome1" ] &&  [ -n "$genome2" ] ; then
       echo "pouet"
       p15=$(sbatch "$DEPENDS"$p12:$p14 "$SCRIPTPATH"/07_dataprep_opt3_opt4_opt5.sh |awk '{print $4}' )
    fi
    opt="synteny_and_Ds"
    echo -e "option for genespace is $opt\n"
    p16=$(sbatch "$DEPENDS"$p15 "$SCRIPTPATH"/08_submit_genespace_paml_and_plot.sh -o "$opt"  |awk '{print $4}' )
fi
################################################################################
# ------section for option 3 : GeneSpace/Synteny + Ds analyses ----------------#
################################################################################
################################################################################
# ----------------- section for option 4 : Ds analyses  -----------------------#
################################################################################
################################################################################
# -------------- section for option 5 : GeneSpace/Synteny analyses ------------#
################################################################################
if [ "$option" == 3 ] || [ "$option" == 4 ] || [ "$option" == 5 ]; then 
       echo -e "option number $option is required \n\n" 
       echo "run data preparation for these option"	
       p15=$(sbatch "$SCRIPTPATH"/07_dataprep_opt3_opt4_opt5.sh |awk '{print $4}' )
fi
if [ "$option" == 3 ] ; then 
    echo "----------------------------------------------------------------"
    echo "       GeneSpace/Synteny + Ds analyses will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------" 
    opt="synteny_and_Ds"
    echo -e "option for genespace is $opt\n"
    p16=$(sbatch "$DEPENDS"$p15 "$SCRIPTPATH"/08_submit_genespace_and_first_plot_only.sh -o "$opt" |awk '{print $4}' )
    p17=$(sbatch "$DEPENDS"$p16 "$SCRIPTPATH"/09_paml_preparation_slurm.sh  |awk '{print $4}' )
    p18=$(sbatch "$DEPENDS"$p17 "$SCRIPTPATH"/10_submit_paml_parallel.sh    |awk '{print $4}' )
    p19=$(sbatch "$DEPENDS"$p18 "$SCRIPTPATH"/11_changepoint_only.sh        |awk '{print $4}' )
    p20=$(sbatch "$DEPENDS"$p19 "$SCRIPTPATH"/12_ideogram_and_circos_after_changepoint.sh  -o "$opt" |awk '{print $4}' )

elif [ "$option" == 4 ] ; then 
    echo "----------------------------------------------------------------"
    echo "         only Ds + associated analyses will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"
    opt="Ds_only"
    p17=$(sbatch "$DEPENDS"$p16 "$SCRIPTPATH"/09_paml_preparation_slurm.sh  |awk '{print $4}' )
    p18=$(sbatch "$DEPENDS"$p17 "$SCRIPTPATH"/10_submit_parallel_paml.sh    |awk '{print $4}' )
    p19=$(sbatch "$DEPENDS"$p18 "$SCRIPTPATH"/11_changepoint_only.sh        |awk '{print $4}' )
    p20=$(sbatch "$DEPENDS"$p19 "$SCRIPTPATH"/12_ideogram_and_circos_after_changepoint.sh  -o "$opt"  |awk '{print $4}' )
elif [ "$option" == 5 ] ; then 
    echo "----------------------------------------------------------------"
    echo "          only GeneSpace/Synteny  analyses will be launched     " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"
    opt="synteny_only"
    p16=$(sbatch "$DEPENDS"$p15 "$SCRIPTPATH"/08_submit_genespace_paml_and_plot.sh -o "$opt"  |awk '{print $4}' )
fi
# Confirm job sumbissions
echo "All jobs submitted"
