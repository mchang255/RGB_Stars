#!/bin/bash
#Preliminary steps (again this is for specific for OH/IR stars project. You may select your own relevant catalogs): Go to the VizieR catalog website. Type "Gaia eDR3" in the search bar and check the boxes for I/350 and I/352 and press "Show table details." Then, check the boxes "I/350/gaiaedr3"  and "I/352/gedr3dis." Then, select "Query selected Tables." Change the search radius from 2 arcmin to 2 arcsec (or some other specified search radius). Check all columns (the main table). On the left hand side, change 50 to unlimited and HTML table to .tsv file (; separated). For the section "Position in," you want to choose the same type as your input coordinates (ex: if you used sexagesimal, select "sexagesimal"). At the top, click "list of targets" where you can upload root-all-final_leThres_coords.csv. After all of that, you can press "submit." When you are finished downloading, give it a name (I chose root-all-final_Thres_eDR3-bailjons.tsv, but again up to you).
#This .sh program takes the output from VizieR .tsv file, fishes out the data for each star and puts them on one line (eDR3 and Bailer-Jones on one line). It then takes stars that have distance values (rpgeo) in the Bailer-Jones catalog and outputs their name, ra, dec, dist(kpc), dist_lower(kpc), dist_higher(kpc), PM, pmRA, e_pmRA, pmDE, e_pmDE, astrometric_excess_noise, astrometric_excess_noise_sig, ipd_gof_harmonic_amplitude, ruwe
#Outputs: ~/mchang-JPL2021/root-all-final_Thres_eDR3-bailjons-hdr.csv - gets the headers from the catalogs selected and puts them on one line
#         ~/mchang-JPL2021/root-all-final_Thres_eDR3-bailjons-cln.csv - The name, original coordinates, and data from .tsv outputted into a .csv (EDR3 and Bailer-Jones data are on one line). Missing data is replaced with -999. If star has multiple data from a catalog, we pair the lines that have the same offset (ex: a star may have two entries of data from EDR3, one data has an offset of 0.4 while the other has an offset of 1.3. There is one entry from Bailer-Jones which has an offset of 0.4. In this case, we put the data that have 0.4 together on one line while the 1.3 data goes onto another line).
#         ~/mchang-JPL2021/root_coords-master.csv - grabs lines that have values for the rgeo (distance) column. Those lines (their names, RA and dec in sexagesimal, rgeo in kpc (along with lower and upper limits), proper motion data, and noise data) will be outputted.
#         ~/mchang-JPL2021/rootx_coords.csv - breaks up the coords-master.csv output into a specified number (x) of files because of ValueErrors returned by get-dict-phot.py when running large files. Right now, I have the program set up so it breaks coords-master.csv into files of 75 lines (so for me that's around 20 files). If you would like to change this, find the line num_lines=75 (near the end of the file) and change it. The outputted files are named from root01 to rootx. These outputs are used for get-dict-phot.py and phot_plotter.py.
#To run this program, if first time running, type chmod +x edr3-bailerjones-dist.sh then enter, then ./edr3-bailerjones-dist.sh root. After that, you can just type ./edr3-bailerjones-dist.sh root
#Ex for running: if I want my root to be 'abc,' then I will type ./edr3-bailerjones-dist.sh abc

root=$1 #assigning root name input to a variable

#creating headers
grep 'SolID' $root\-all-final_Thres_eDR3-bailjons.tsv | tail -1 > $root\-all-final_Thres_eDR3-bailjons-hdr.csv

bail=$(grep 'rgeo' $root\-all-final_Thres_eDR3-bailjons.tsv | tail -1)

sed -i '' -e "s/\$/;$bail/g" $root\-all-final_Thres_eDR3-bailjons-hdr.csv

#how many columns there are for each catalog
length_edr3=$(grep 'SolID' $root\-all-final_Thres_eDR3-bailjons.tsv | tail -1 | awk 'BEGIN{FS=";"}{print NF}')

length_bail=$(grep 'rgeo' $root\-all-final_Thres_eDR3-bailjons.tsv | tail -1 | awk 'BEGIN{FS=";"}{print NF}')

header=$(cat $root\-all-final_Thres_eDR3-bailjons-hdr.csv)

units_gaia=$(grep 'deg' $root\-all-final_Thres_eDR3-bailjons.tsv | tail -2 | head -1 | sed 's/ ;/--;/g')
units_bail=$(grep 'deg' $root\-all-final_Thres_eDR3-bailjons.tsv | tail -1 | sed 's/ ;/--;/g')

#ouputting header and units to top of cln file
echo "Name;Input RA;Input Dec;$header" > $root\-all-final_Thres_eDR3-bailjons-cln.csv
echo "--;\"h:m:s\";\"h:m:s\";$units_gaia;$units_bail" >> $root\-all-final_Thres_eDR3-bailjons-cln.csv

n=$(wc -l < $root\-all-final_leThres_coords.csv)

#Section 1: This section extracts data from the .tsv into the cln.csv
for i in $(seq 1 $n)
do
    #Section 1a: In the .tsv file, for each source, it will say Constraint: input RA input dec. This section of the code extracts data that goes from one Constraint to the next Constraint. However, there are some sources that have data for one catalog but not for another. Some sources might have more data in one catalog than the other. Therefore, we need to account for this.
    let s=$i+1
    a=$(head -$i $root\-all-final_leThres_coords.csv | tail -1)
    b=$(head -$s $root\-all-final_leThres_coords.csv | tail -1)
    name=$(awk -v i="$s" 'BEGIN{FS=", ";OFS=", "} NR==i {print $1}' $root\-all-final_leThres.csv)
    ra=$(awk -v i="$s" 'BEGIN{FS=", ";OFS=", "} NR==i {print $4}' $root\-all-final_leThres.csv)
    dec=$(awk -v i="$s" 'BEGIN{FS=", ";OFS=", "} NR==i {print $5}' $root\-all-final_leThres.csv)
    
    gaia=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/350\/gaiaedr3/,/I\/352\/gedr3dis/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' )
    bailer=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/352\/gedr3dis/,/#Constraint $b/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g')
    number1=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/350\/gaiaedr3/,/I\/352\/gedr3dis/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | wc -l)
    number2=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/352\/gedr3dis/,/#Constraint $b/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | wc -l)
    #Section 1a finished.
    
    #Section 1b: Let's say that for one source it has no data for one catalog. In this case, for the catalog that does not have data, we replace the blanks with -999. For the catalog that does have data, we need to see how many entries there are. If it has just one entry, then we can immediately print the data onto one line. However, if there are multiple entries, we need to print multiple lines, one of for each entry.
    if [[ $number1 -eq 0 || $number2 -eq 0 ]]
    then
        if [[ $gaia == "" ]]
        then
            gaia=$(python -c "print('-999;' * $length_edr3)" | sed 's/.$//')
        fi

        if [[ $bailer == "" ]]
        then
            bailer=$(python -c "print('-999;' * $length_bail)" | sed 's/.$//')
        fi

        if [ $number1 -gt 1 ]
        then
            for u in $(seq 1 $number1)
            do
                gaia=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/350\/gaiaedr3/,/I\/352\/gedr3dis/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | head -$u | tail -1)
                echo "$name;$ra;$dec;$gaia;$bailer"
            done
        elif [ $number2 -gt 1 ] #replaced $number2 -gt 2 with $number2 -gt 1
        then
            for u in $(seq 1 $number2)
            do
                bailer=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/352\/gedr3dis/,/#Constraint $b/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | head -$u | tail -1 )
                echo "$name;$ra;$dec;$gaia;$bailer"
            done
        else
            echo "$name;$ra;$dec;$gaia;$bailer"
        fi
    #Section 1b finished.
    else
    #Section 1c: Let's say now that one catalog has more entries than another. In this case, we look for the catalog that has the greater number of entries. We look at the different offsets (the first column of the data). Then, for the all the data in the source's section, we grab the lines (data) that have the same offset. If the number of lines that have the same offset is even (meaning the value of the offset is present in both gaia and bailer-jones), we will put these together on one line. If the number of lines is odd (really just one line), we put this data on another line along with -999 for the other catalog as that offset does not exist for the other catalog. If the two catalogs have the same entries as each other, we put one gaia and one bailer-jones on one line.
        if [ $number1 -gt $number2 ]
        then
            for u in $(seq 1 $number1)
            do
                offset=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/350\/gaiaedr3/,/I\/352\/gedr3dis/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | head -$u | tail -1 | awk 'BEGIN{FS=";";OFS=";"}{print $1}')
                is_even=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}' | sed 's/ ;/-999;/g' | grep "^$offset" | wc -l)
                if [ $(expr $is_even % 2) -eq 0 ]
                then
                    offset_lines=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}' | sed 's/ ;/-999;/g' | grep "^$offset" | tr '\n' ';')
                else
                    filler=$(python -c "print('-999;' * $length_bail)" | sed 's/.$//')
                    offset_lines="$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/350\/gaiaedr3/,/I\/352\/gedr3dis/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}' | sed 's/ ;/-999;/g' | grep "^$offset" | tr '\n' ';')$filler"
                fi
                echo "$name;$ra;$dec;$offset_lines"
            done
        elif [ $number2 -gt $number1 ]
        then
            for u in $(seq 1 $number2)
            do
                offset=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/352\/gedr3dis/,/#Constraint $b/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | head -$u | tail -1 | awk 'BEGIN{FS=";";OFS=";"}{print $1}')
                is_even=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}' | sed 's/ ;/-999;/g' | grep "^$offset" | wc -l)
                if [ $(expr $is_even % 2) -eq 0 ]
                then
                    offset_lines=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}' | sed 's/ ;/-999;/g' | grep "^$offset" | tr '\n' ';')
                else
                    filler=$(python -c "print('-999;' * $length_edr3)")
                    offset_lines="$filler$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/352\/gedr3dis/,/#Constraint $b/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}' | sed 's/ ;/-999;/g' | grep "^$offset" | tr '\n' ';')"
                fi
                echo "$name;$ra;$dec;$offset_lines"
            done 
        elif [ $number1 -eq $number2 ]
        then
            for u in $(seq 1 $number1)
            do
                gaia_line=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/350\/gaiaedr3/,/I\/352\/gedr3dis/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | head -$u | tail -1)
                bailer_line=$(sed -n "/#Constraint $a/,/#Constraint $b/p" $root\-all-final_Thres_eDR3-bailjons.tsv | sed -n "/I\/352\/gedr3dis/,/#Constraint $b/p" | awk 'BEGIN{FS=";"}/^[0-9]/{print $0}'| sed 's/ ;/-999;/g' | head -$u | tail -1)
                echo "$name;$ra;$dec;$gaia_line;$bailer_line"
            done            
        fi
    fi
done >> $root\-all-final_Thres_eDR3-bailjons-cln.csv
#Section 1c finished.
#Section 1 finished.

sed -i '' -e 's/;/,/g' $root\-all-final_Thres_eDR3-bailjons-cln.csv

#extracting sources that have an rgeo value. creating the coords file which has the name, RA and dec (sexagesimal), distance in kpc (upper and lower limits), PM, pnRA, e_pmRA, pmDE, e_pmDE, astrometric_excess_noise, astrometric_excess_noise_sig, ipd_gof_harmonic_amplitude, ruwe
sed '1,2d' $root\-all-final_Thres_eDR3-bailjons-cln.csv | awk 'BEGIN{FS=",";OFS=","} {if($128!=-999)print $1,$5,$6,($128 / 1000),($129 / 1000),($130 / 1000),$19,$20,$21,$22,$23,$40,$41,$58,$62}' > $root\_coords-master.csv

#breaking the file up into x number of files
num_lines=$(wc -l < $root\_coords-master.csv)
#num_files=$(echo "( ${num_lines} + 9 ) / 10" | bc)
num_lines=75
split -l "$num_lines" $root\_coords-master.csv temp_name

o=1
for file in temp_name*
do
 if [ $o -le 9 ]
 then
     mv $file $root\0${o}_coords.csv
 else
     mv $file $root${o}_coords.csv
 fi
 o=$(( o + 1 ))
done
