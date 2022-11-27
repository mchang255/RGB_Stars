#!/bin/bash
#Preliminary steps (Note this is for OH/IR stars, so if there are other star types you want, YOU need to come up with your own steps):  
#Step 1: Type "OH/IR stars" into VizieR (https://vizier.u-strasbg.fr/viz-bin/VizieR):
#           -When viewing the catalogs, some catalogs will be missing information. Also, some catalogs will contain subparts called tables. In order for the catalog to be used, or part (tables) of the catalog, it NEEDS to contain name, right ascension, and declination. To check if the catalog satisfies that condition, click the hyperlink (it starts with J/..). On that page, there will be a submit button, right next to the reset all button. Click the submit button. You should see the name, right ascension, and declination within the first several columns. If one of those is missing for one table, but is visible in other tables, the catalog is good. If all tables/catalog are missing at least one of those features, you may have to do more digging and view the "Read Me" file on the query page to go further and find the original catalog. Otherwise, DO NOT use that catalog. The following catalogs apply (aka DO NOT USE):
#           -Amiri+, 2012 (J/A+A/538/A136)
#           -Lewis 1997 (J/ApJS/109/489). However, when looking at the Read Me, the stars are actually taken from the Arecibo 1612 mission, so enter "Arecibo 1612" in the search bar on the main page and query those (Lewis 1994 J/ApJS/93/549).
#           -Lepine+, 1995 (J/A+A/299/453)
#           Also, don't use the following catalogs because they are included in others:
#           -Chengalur+ 1993 (J/ApJS/89/189) (included in Engels+ 2015)
#           -Chen+, 2001 (J/A+A/368/1006) (included in Engels+ 2015)
#           Next, when it comes to right ascension and declination, we need to have the coordinates in both degrees and sexagesimal. With each catalog, and even each table, the coordinates are given in one of those two. The other coordinates need to be supplied by VizieR. So, with each catalog, check if the coordinates are given in sexagesimal (h:m:s) or degrees. After checking, if the coordinates by the catalog are given in sexagesimal, then we need to select decimal for VizieR to supply (for more on this, we will see in the next following steps). Sometimes, the tables within each catalog have coordinates different from each other (meaning some tables will be given in sexagesimal and some in degrees, all within in the same catalog). In that case, you will have to go table by table and download the information. Now, we move on to downloading the information. Download everything (ex: photometry, coordinates) about star, so check all columns before pressing submit. Change 50 to unlimited and HTML table to .tsv file (tab-separated) on the left hand side. For where it says "Position in:", choose the coordinate type that IS NOT given by the catalog. Repeat for remaining catalogs. Rename and place all of these .tsv files into a folder (you can name it whatever you want, but the name of the folder will be the root name for all of the outputs from running search.sh) (location ~/mchang-JPL2021/root). 
#           To make things clearer, I will give an example of the above. Let me take the Sevenster catalog. When I see the information displayed about Sevenster, I see that for the first table (J/AJ/123/2772/AGB), the coordinates given are in degrees. For the other two tables, (J/AJ/123/2772/doublede and J/AJ/123/2772/tablevla), the coordinates are given in sexagesimal. So, because of this, I will have to split up the catalog into tables that have degrees and tables that have sexagesimal. So, J/AJ/123/2772/AGB will be one file and J/AJ/123/2772/doublede and J/AJ/123/2772/tablevla will be in another file. For J/AJ/123/2772/AGB, when I'm downloading this file, I will select "sexagesimal" (where it says "Position in:", more on this below). Likewise J/AJ/123/2772/doublede and J/AJ/123/2772/tablevla, I will select "decimal."
#           So far, the catalogs I have are:
#           -J/ApJS/93/549
#           -J/AJ/123/2772
#           -J/A+AS/92/43
#           -J/A+AS/127/185
#           -J/A+A/391/967
#           -J/A+A/446/773
#           -J/AJ/127/501
#           -J/A+AS/128/35
#           -J/A+A/431/779
#           -J/A+A/562/A68
#               -by extension of above, J/A+A/562/A68/maps
#           -J/MNRAS/505/6051
#This code runs through all the files in the 'root' folder, which contains all the various catalogs of root stars downloaded from VizieR, and extracts the name of the star, RA (both deg and sexag), dec (both deg and sexag), name of the catalog, and code of the catalog from each file. A new file is created for each catalog and then stored in the subdirectory 'root-clean'. Then, we concatenate all the files and query the stars in Simbad to get even more accurate coordinates, alternate names, and the differences between the input and output coordinates. Then, we select stars that are below a certain offset. See below for outputs and what differentiate them.
#Outputs: ~/mchang-JPL2021/root-all-raw1.csv - master list, not sorted
#               ~/mchang-JPL2021/root-all-raw2.csv - master list sorted by RA
#               ~/mchang-JPL2021/root-all-raw3.csv - master list (identical duplicates of name, RA, and dec are deleted from raw2), star coordinates are entered into Simbad, better coordinates (both decimal and sexagesimal) are returned, non-IRAS named stars replaced with IRAS names(if possible), alternate names are sought, the search radius of Simbad, and the offset between input and output coordinates are calculated. What this code does: the coordinates are entered into Simbad starting at a search radius of 2 arcseconds. If an object is found, the name of the object (IRAS name, if not, the main name listed on Simbad in bold), the coordinates (sexagesimal and deg) from Simbad, alternate names, the offset, and search radius are printed. If no object is found, we increase  the search radius by 1 arcsec until an object is found.
#               ~/mchang-JPL2021/root-all-final.csv - After root-all-raw3.csv, we will delete stars that are identical duplicates in name, RA (deg), and dec (deg)
#               ~/mchang-JPL2021/root-all-final_rejects.csv - The stars that were deleted from root-all-raw3.csv will be outputted to this file
#               ~/mchang-JPL2021/root-all-final_leTHRES.csv - Extracting stars that were found less than number arcseconds from the given coordinates.
#               ~/mchang-JPL2021/root-all-final_leTHRES_rejects.csv - The stars that were not chosen from root-all-final.csv will be outputted to this file
#               ~/mchang-JPL2021/root-all-final_leTHRES_coords.csv - Extracting just the columns that have the new RA and dec in sexagesimal coordinates (this can be inputted in VizieR)
#Every output acts on the previous (ex: root-all-raw2.csv acts on root-all-raw1.csv, root-all-raw3.csv acts on root-all-raw2.csv)
#To run this code, type ./search.sh root threshold into terminal (if it's the first time or if just typing ./search.sh root threshold gives a permission denied error, type chmod +x search.sh, hit enter, then ./search.sh root threshold). Threshold is in arcsec. If you are worried about your computer falling asleep while the code is running, then type: caffeinate ./search.sh root threshold
#Ex of running: I want my root name to be 'abc', and my threshold to be 10 arcseconds. In the terminal, I will type ./search.sh abc 10 (or caffeinate ./search.sh abc 10).

root=$1
thresh=$2

mkdir ~/mchang-JPL2021/$root/$root\-clean

#Section 1. This section of the code runs through all the files in the 'root' folder, which contains all the various catalogs of OH-IR stars, and extracts the name of the star, RA, dec, name of the catalog, and code of the catalog from each file. A new file is created for each catalog and then stored in the folder 'root-clean'. Still not done adding catalogs to list and if there are more in the future, you are going to have to add more code.
for file in ~/mchang-JPL2021/$root/*
do 
    if [[ $(basename "$file") == "arecibo1612-1.tsv" ]]
    then
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$4","$1","$2","$26","$27}' $file | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/arecibo1612-final.csv
    elif [[ $(basename "$file") == "arecibo1612-2.tsv" ]]
    then
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$4","$9","$10","$1","$2}' $file | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' >> ~/mchang-JPL2021/$root/$root\-clean/arecibo1612-final.csv
    elif [[ $(basename "$file") == "jim-est05.tsv" ]]
    then
        awk 'BEGIN{FS="\t";OFS="\t"} NR < 437 {print $0}' $file | awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$3","$1","$2","$4","$5}' | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/jim-est05-final.csv
    elif [[ $(basename "$file") == "jim-est06.tsv" ]]
    then
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$4","$1","$2","$6","$7}' $file | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/jim-est06-final.csv
    elif [[ $(basename "$file") == "lewis04.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$3","$16","$17","$1","$2}' $file | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/lewis04-final.csv
    elif [[ $(basename "$file") == "lindqvist.tsv" ]]
    then
        awk 'BEGIN{FS="\t"} /^[0-9]/{print $4","$1","$2","$20","$21}' $file | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/lindqvist-final.csv
    elif [[ $(basename "$file") == "nyman.tsv" ]]
    then
        awk 'BEGIN{FS="\t";OFS="\t"} NR < 211 {print $0}' $file | awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$3","$1","$2","$13","$14}' | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/nyman-final.csv
    elif [[ $(basename "$file") == "sevenster1.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$7","$30","$31","$1","$2}' $file | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/sevenster-final.csv
    elif [[ $(basename "$file") == "sevenster2.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$10","$1","$2","$4","$5}' $file | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' >> ~/mchang-JPL2021/$root/$root\-clean/sevenster-final.csv
    elif [[ $(basename "$file") == "sevenster3.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$18","$1","$2","$6","$7}' $file | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' >> ~/mchang-JPL2021/$root/$root\-clean/sevenster-final.csv
    elif [[ $(basename "$file") == "sjouwerman02-1.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "OH"$3","$1","$2","$6","$7}' $file | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/sjouwerman02-final.csv
    elif [[ $(basename "$file") == "sjouwerman02-2.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "OH"$3","$1","$2","$5","$6}' $file | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' >> ~/mchang-JPL2021/$root/$root\-clean/sjouwerman02-final.csv
    elif [[ $(basename "$file") == "sjouwerman98.tsv" ]]
    then
        awk 'BEGIN{FS="\t"} /^[0-9]/{print "OH"$4","$1","$2","$5","$6}' $file | awk -v title="$(grep -m1 '#Title:' $file)" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/sjouwerman98-final.csv
    elif [[ $(basename "$file") == "engels1.tsv" ]]
    then
       title="$(grep -m1 '#Title:' $file)"
       title="${title//,/}"
       awk 'BEGIN{FS="\t"} !/</{print $0}' $file | awk 'BEGIN{FS="\t"}/^[0-9]/{print $7","$1","$2","$10","$11}' | sed 's/IRAS/IRAS /g' | sed 's/IRAS  /IRAS /g' | awk -v title="$title" -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/engels1-final.csv
    elif [[ $(basename "$file") == "engels2.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t"} !/</{print $0}' $file | awk 'BEGIN{FS="\t"}/^[0-9]/{print $6","$1","$2","$7","$8}' | awk -v title="$title" -v name="$(grep -m1 '#Table' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,title,name}' | sed 's/#Title: //g;s/#Table//g' > ~/mchang-JPL2021/$root/$root\-clean/engels2-final.csv
    elif [[ $(basename "$file") == "jim-est21.tsv" ]]
    then
        title="$(grep -m1 '#Title:' $file)"
        title="${title//,/}"
        awk 'BEGIN{FS="\t";OFS="\t"} NR < 496 {print $0}' $file | awk 'BEGIN{FS="\t"} /^[0-9]/{print "IRAS",$3","$5","$6","$1","$2}' | awk -v name="$(grep -m1 '#Name:' $file)" 'BEGIN{FS=",";OFS=","}{print $1,$2,$3,$4,$5,"Arecibo sample variability properties (Jimenez-Esteban+ 2021)",name}' | sed 's/#Title: //g;s/#Name: //g' > ~/mchang-JPL2021/$root/$root\-clean/jim-est21-final.csv
    fi
done
#Section 1 finished.

#Section 2: This section of the code concatenates all the files in 'root-clean' to create one master list called 'root-all-raw1.csv' (no duplicates have been deleted yet). It also gets rid of alignment issues in the first column. It also makes the objects names all uppercase and changes the field separator from comma to comma + space. A header is added to $root\-all-raw1.csv. Then, root-all-raw1.csv sorted in ascending RA order to make root-all-raw2.csv. Then identical duplicates in name, RA (deg), and Dec (deg) are deleted, extra white spaces in the second and third column are deleted to make root-all-raw3.csv. A different header is then used for root-all-raw3.csv.
cat ~/mchang-JPL2021/$root/$root\-clean/* | awk 'BEGIN{FS=",";OFS=","}{gsub(/[ \t]+$/, "", $1); print $0}' | awk 'BEGIN{FS=",";OFS=","} $1 = toupper($1)' | awk 'BEGIN{FS=",";OFS=", "} {$1=$1}1' > $root\-all-raw1.csv

sed -i '' -e '1 i\
#Name, RA (deg), Dec (deg), RA (sexag), Dec (sexag), Vizier Cat, Ref Number' $root\-all-raw1.csv

sort -nk2 -t"," $root\-all-raw1.csv > $root\-all-raw2.csv

cat $root\-all-raw2.csv | awk 'BEGIN{FS=", ";OFS=", "} !seen[$1,$2,$3]++' | awk 'BEGIN{FS=", ";OFS=", "}{gsub(" ", "", $2); print $0}' | awk 'BEGIN{FS=", ";OFS=", "}{gsub(" ", "", $3); print $0}' > $root\-all-raw3.csv

sed -i '' -e '1s/.*/#Name, New RA (deg), New Dec (deg), New RA (sexag), New Dec (sexag), Vizier Cat, Ref Number, Alt Names, Orig RA (sexag), Orig Dec (sexag), Offset (arcsec), Search Radius (arcsec), Original Name, Comment on Name/g' $root\-all-raw3.csv
#Section 2 finished.

#Section 3: This while loop loops line by line of root-all-raw3.csv. and goes to Simbad to fetch a variety of things using the coordinates given from root-all-raw3.csv. The components of the loop will be broken into subsections.
n=$(wc -l < $root\-all-raw3.csv)
i=2

while [ $i -le $n ]
do
    #Section 3a: Line by line, this section extracts the name, RA and dec (both deg and sexag). Then, the link to the Simbad is provided using the RA and dec (deg). As for what the weblink looks like, it depends whether the declination has a + or - and that is because how the web browser formats +, -, and blank spaces (a blank space is + in the address bar while + is +%2B, so it's important we get the correct formatting). The weblink starts with the search radius 2 arcseconds.
    first_name=$(awk 'BEGIN{FS=", "}{print $1}' $root\-all-raw3.csv | head -$i | tail -1)
    ra=$(awk 'BEGIN{FS=", "}{print $2}' $root\-all-raw3.csv | head -$i | tail -1)
    dec=$(awk 'BEGIN{FS=", "}{print $3}' $root\-all-raw3.csv | head -$i | tail -1)
    ra_s=$(awk 'BEGIN{FS=", "}{print $4}' $root\-all-raw3.csv | head -$i | tail -1)
    dec_s=$(awk 'BEGIN{FS=", "}{print $5}' $root\-all-raw3.csv | head -$i | tail -1)
    
    if [[ $dec == *+* ]]
    then
        new_dec="${dec//+/+%2B}"
        name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=")
        echo "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList="
        loop="https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList="
    else
        name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=")
        echo "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList="
        loop="https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList="
    fi
    #Section 3a finished.
    
    #Section 3b: On the first try (starting with 2 arcseconds), we may get multiple sources. If that's the case, we look at the following section. The section of the code then loops through each source and checks if the source has an IRAS name. The first entry to have an IRAS name is returned. If none of the sources have an IRAS name, then we use the source that has the least offset, and we choose the first name that is provided by Simbad (if you search a source via coordinates, we take the big bolded name at the top of the webpage). From there -- and again because of formatting, + and - declinations take separate paths -- , we get the alternate names and the Simbad-provided RA and dec (both decimal and sexag). Sometimes, the RA and dec are not given in decimal, in that case we need to take the RA and dec in sexagesimal and go to a converter (I chose http://www.astrouw.edu.pl/~jskowron/ra-dec/ as my converter). Then, sometimes the code might behave weirdly in that the final RA and dec may be blank or contain a value that is not expected, in which case, if they do, the input source is looped again to get a different, more correct value. This serves as a check. Assuming the value is correct, we then replace the input name with the Simbad-provided name, same with RA and dec. We place the RA and dec in sexag and the original in other columns, so we can compare. We also note what the offset between the source from Simbad and the input source and also the search radius (in this case, it's 2 arcseconds). Then, we compare the original and Simbad-provided names. If the Simbad name is IRAS and the original is not, we say that an IRAS name is found. If the Simbad name is an IRAS name, and it matches the original, we say No issues. If the Simbad name is an IRAS name, and the original is also an IRAS name, but the two are different, we say that the original is incorrect, and that we have replaced it. If the Simbad name is not an IRAS name, we say that no IRAS name has been found.
    if [[ $name == *Number\ of\ rows* ]]
    then
        rows=$(curl --silent "${loop}" | grep 'simbad.u-strasbg.fr/simbad/sim-id' | grep 'submit=submit' | sed 's/.*HREF="\/\///;s/".*//')
        row_number=$(curl --silent "${loop}" | grep 'simbad.u-strasbg.fr/simbad/sim-id' | grep 'submit=submit' | sed 's/.*HREF="\/\///;s/".*//' | wc -l)
        start_row=1

        while [ $start_row -le $row_number ]
        do
            link=$(echo "${rows}" | head -$start_row | tail -1)
            webpage=$(curl --silent "${link}")
            echo $link
            if [[ $webpage == *IRAS* ]]
            then
                replace_name=$(curl --silent "${link}" | grep "IRAS [A-Z0-9]" | sed 's/.*"//;s/<.*//;s/>//' | grep "IRAS [A-Z0-9]" | tail -1)
                altname=$(curl --silent "${link}" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                
                altname="${altname//,/:}"
                if [[ $dec == \+* ]]
                then
                    new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
                    dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
                    new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
                    new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)
                    
                    if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                    then
                        dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                        new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                        dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                        new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                        new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                    fi
                else
                    new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
                    dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
                    new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
                    new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1)
                    
                    if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                    then
                        dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                        new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                        dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                        new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                        new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
                    fi
                fi
                
                if [[ $new_ra == "271.087458" ]] || [[ $dec3 == "-29.519319" ]]
                then
                    i=$(( i-1 ))
                elif [[ $replace_name == "" ]] || [[ $new_ra == "" ]] || [[ $dec3 == "" ]] || [[ $new_ra_s == "" ]] || [[ $new_dec_s == "" ]] || [[ $altname == "" ]]
                then
                    i=$(( i-1 ))
                else
                    
                    offset=$(python offset.py $ra $dec $new_ra $dec3)
                    offset="${offset//arcsec/}"
                    
                    sed -i '' -e "${i}s/$first_name/$replace_name/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$ra/$new_ra/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$dec/$dec3/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$ra_s/$new_ra_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$dec_s/$new_dec_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $altname/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $ra_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $dec_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $offset/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, 2/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $first_name/g" $root\-all-raw3.csv
                    echo "Several rows: replace ${first_name} with ${replace_name}, ${ra} with ${new_ra}, ${dec} with ${dec3}"

                    if [[ $replace_name == *IRAS* ]]
                    then
                        if [[ $first_name !=  *IRAS* ]]
                        then
                            sed -i '' -e "${i}s/\$/, ${first_name} replaced with IRAS in Col 1./" $root\-all-raw3.csv
                        else
                            if [[ $first_name == $replace_name ]]
                            then
                                sed -i '' -e "${i}s/\$/, No issues./" $root\-all-raw3.csv
                            else
                                sed -i '' -e "${i}s/\$/, ${first_name} was INCORRECT. Corrected in Col 1./" $root\-all-raw3.csv
                            fi
                        fi
                    else
                        sed -i '' -e "${i}s/\$/, Coord search did not find IRAS ID./" $root\-all-raw3.csv
                    fi
                fi

                break
            else
                if [ $start_row -eq $row_number ]
                then
                    link=$(echo "${rows}" | head -1)
                    echo $link
                    replace_name=$(curl --silent "${link}" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//;')
                    altname=$(curl --silent "${link}" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')

                    if [[ $dec == \+* ]]
                    then
                        new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
                        dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
                        new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
                        new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)  
                        
                        if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                        then
                            dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                            new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                            dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                            new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                            new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                        fi
                    else
                        new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
                        dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
                        new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
                        new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1)
                        
                        if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                        then
                            dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                            new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                            dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                            new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                            new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
                        fi
                    fi
                    
                    if [[ $new_ra == "271.087458" ]] || [[ $dec3 == "-29.519319" ]]
                    then
                        i=$(( i-1 ))
                    elif [[ $replace_name == "" ]] || [[ $new_ra == "" ]] || [[ $dec3 == "" ]] || [[ $new_ra_s == "" ]] || [[ $new_dec_s == "" ]] || [[ $altname == "" ]]
                    then
                        i=$(( i-1 ))
                    else
                        offset=$(python offset.py $ra $dec $new_ra $dec3)
                        offset="${offset//arcsec/}"
                        
                        sed -i '' -e "${i}s/$first_name/$replace_name/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/$ra/$new_ra/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/$dec/$dec3/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/$ra_s/$new_ra_s/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/$dec_s/$new_dec_s/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/\$/, $altname/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/\$/, $ra_s/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/\$/, $dec_s/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/\$/, $offset/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/\$/, 2/g" $root\-all-raw3.csv
                        sed -i '' -e "${i}s/\$/, $first_name/g" $root\-all-raw3.csv
                        echo "Several rows: replace ${first_name} with ${replace_name}, ${ra} with ${new_ra}, ${dec} with ${dec3}"
                        
                        if [[ $replace_name == *IRAS* ]]
                        then
                            if [[ $first_name !=  *IRAS* ]]
                            then
                                sed -i '' -e "${i}s/\$/, ${first_name} replaced with IRAS in Col 1./" $root\-all-raw3.csv
                            else
                                if [[ $first_name == $replace_name ]]
                                then
                                    sed -i '' -e "${i}s/\$/, No issues./" $root\-all-raw3.csv
                                else
                                    echo "Andrei knows everything"
                                    sed -i '' -e "${i}s/\$/, ${first_name} was INCORRECT. Corrected in Col 1./" $root\-all-raw3.csv
                                fi
                            fi
                        else
                            sed -i '' -e "${i}s/\$/, Coord search did not find IRAS ID./" $root\-all-raw3.csv
                        fi
                    fi
                    break
                else
                    start_row=$(( start_row+1 ))
                fi
            fi
        done
    #Section 3b finished
    
    #Section 3c: On the first try (starting with 2 arcseconds), we may get "No astronomical is found." If that's the case, we look at the following section. The section of the code then keeps increasing the search radius by 1 arcsecond until a source is found. However, we may get the situation described in section 3b, which is upon our first return, there may be multiple entries. In which case, I copied and pasted the code from section 3b. If just one source is returned, the section of the code breaks into whether the source has an IRAS name or not. From there -- and again because of formatting, + and - declinations take separate paths -- , we get the alternate names and the Simbad-provided RA and dec (both decimal and sexag). Sometimes, the RA and dec are not given in decimal, in that case we need to take the RA and dec in sexagesimal and go to a converter (I chose http://www.astrouw.edu.pl/~jskowron/ra-dec/ as my converter). Then, sometimes the code might behave weirdly in that the final RA and dec may be blank or contain a value that is not expected, in which case, if they do, the input source is looped again to get a different, more correct value. This serves as a check. Assuming the value is correct, we then replace the input name with the Simbad-provided name, same with RA and dec. We place the RA and dec in sexag and the original in other columns, so we can compare. We also note what the offset between the source from Simbad and the input source and also the search radius at which the object was returned. Then, we compare the original and Simbad-provided names. If the Simbad name is IRAS and the original is not, we say that an IRAS name is found. If the Simbad name is an IRAS name, and it matches the original, we say No issues. If the Simbad name is an IRAS name, and the original is also an IRAS name, but the two are different, we say that the original is incorrect, and that we have replaced it. If the Simbad name is not an IRAS name, we say that no IRAS name has been found.    
    elif [[ $name == *No\ astronomical* ]]
    then
#         arcsec=6
        start_sec=3
        while true
        do
            if [[ $dec == \+* ]]
            then
                second_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=")
                echo "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList="
                loop="https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList="
            else
                second_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=")
                echo "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList="
                loop="https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList="
            fi
            
            if [[ $second_name == *No\ astronomical* ]]
            then
#                 if [ $start_sec -eq $arcsec ]
#                 then
#                     sed -i '' -e "${i}s/\$/, , No crossmatch found in Simbad\/misID'd IRAS source/" $root\-all-raw3.csv
#                     echo "${first_name}, ${ra} not found in Simbad"
#                 fi
                start_sec=$(( start_sec+1 ))
            else
                if [[ $second_name == *Number\ of\ rows* ]]
                then
                    rows=$(curl --silent "${loop}" | grep 'simbad.u-strasbg.fr/simbad/sim-id' | grep 'submit=submit' | sed 's/.*HREF="\/\///;s/".*//')
                    row_number=$(curl --silent "${loop}" | grep 'simbad.u-strasbg.fr/simbad/sim-id' | grep 'submit=submit' | sed 's/.*HREF="\/\///;s/".*//' | wc -l)
                    start_row=1
                        
                    while [ $start_row -le $row_number ]
                    do
                        link=$(echo "${rows}" | head -$start_row | tail -1)
                        echo $link
                        webpage=$(curl --silent "${link}")
                        if [[ $webpage == *IRAS* ]]
                        then
                            replace_name=$(curl --silent "${link}" | grep "IRAS [A-Z0-9]" | sed 's/.*"//;s/<.*//;s/>//' | grep "IRAS [A-Z0-9]" | tail -1)
                            altname=$(curl --silent "${link}" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                            altname="${altname//,/:}"
                            if [[ $dec == \+* ]]
                            then
                                new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
                                dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
                                new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
                                new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)  
                                
                                if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                                then
                                    dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                                    new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                                    dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                                    new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                                    new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                                fi
                            else
                                new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
                                dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
                                new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
                                new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1)  
                                
                                if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                                then
                                    dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                                    new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                                    dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                                    new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                                    new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
                                fi
                            fi
                                
                            break
                        else
                            if [ $start_row -eq $row_number ]
                            then
                                link=$(echo "${rows}" | head -1)
                                echo $link
                                replace_name=$(curl --silent "${link}" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//;')
                                altname=$(curl --silent "${link}" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                                altname="${altname//,/:}"
                                    
                                if [[ $dec == \+* ]]
                                then
                                    new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
                                    dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
                                    new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
                                    new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)  
                                    
                                    if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                                    then
                                        dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                                        new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                                        dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                                        new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                                        new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                                    fi
                                else
                                    new_ra=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
                                    dec3=$(curl --silent "${link}" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
                                    new_ra_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
                                    new_dec_s=$(curl --silent "${link}" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1)  
                                    
                                    if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                                    then
                                        dec_ra=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                                        new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                                        dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                                        new_ra_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                                        new_dec_s=$(curl --silent "${link}" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
                                    fi
                                fi
                                break
                            else
                                start_row=$(( start_row+1 ))
                            fi
                        fi
                    done
                elif [[ $second_name == *IRAS* ]]
                then
                    if [[ $dec == \+* ]]
                    then
                        replace_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep  "IRAS [A-Z0-9]" | sed 's/.*"//;s/<.*//;s/>//' | grep  "IRAS [A-Z0-9]" | tail -1)
                        new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
                        dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
                        new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
                        new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)  
                        altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                        altname="${altname//,/:}"
                        
                        if [[ $replace_name != IRAS\ * ]]
                        then
                            replace_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//')
                        fi
                        
                        if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                        then
                            dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                            new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                            dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                        fi
                    else
                        replace_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep  "IRAS [A-Z0-9]" | sed 's/.*"//;s/<.*//;s/>//' | grep  "IRAS [A-Z0-9]" | tail -1)
                        new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
                        dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
                        new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
                        new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1) 
                        altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                        altname="${altname//,/:}"
                        
                        if [[ $replace_name != IRAS\ * ]]
                        then
                            replace_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//')
                        fi
                        
                        if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                        then
                            dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                            new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                            dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
                        fi
                    fi
                else
                    if [[ $dec == \+* ]]
                    then
                        replace_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//;')
                        new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
                        dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
                        new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
                        new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)  
                        altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                        altname="${altname//,/:}"
                        
                        if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                        then
                            dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                            new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                            dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                        fi
                    else
                        replace_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//;')
                        new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
                        dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
                        new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
                        new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1)  
                        altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
                        altname="${altname//,/:}"
                        
                        if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
                        then
                            dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                            new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                            dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=${start_sec}&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
                        fi
                    fi
                fi
                
                if [[ $new_ra == "271.087458" ]] || [[ $dec3 == "-29.519319" ]]
                then
                    i=$(( i-1 ))
                elif [[ $replace_name == "" ]] || [[ $new_ra == "" ]] || [[ $dec3 == "" ]] || [[ $new_ra_s == "" ]] || [[ $new_dec_s == "" ]] || [[ $altname == "" ]]
                then
                    i=$(( i-1 ))
                else
                    offset=$(python offset.py $ra $dec $new_ra $dec3)
                    offset="${offset//arcsec/}"
                    
                    sed -i '' -e "${i}s/$first_name/$replace_name/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$ra/$new_ra/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$dec/$dec3/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$ra_s/$new_ra_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/$dec_s/$new_dec_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $altname/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $ra_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $dec_s/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $offset/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $start_sec/g" $root\-all-raw3.csv
                    sed -i '' -e "${i}s/\$/, $first_name/g" $root\-all-raw3.csv
                    echo "No astronomical found initially: replace ${first_name} with ${replace_name}, ${ra} with ${new_ra}, ${dec} with ${dec3}"

                    if [[ $replace_name == *IRAS* ]]
                    then
                        if [[ $first_name !=  *IRAS* ]]
                        then
                            sed -i '' -e "${i}s/\$/, ${first_name} replaced with IRAS in Col 1./" $root\-all-raw3.csv
                        else
                            if [[ $first_name == $replace_name ]]
                            then
                                sed -i '' -e "${i}s/\$/, No issues./" $root\-all-raw3.csv
                            else
                                sed -i '' -e "${i}s/\$/, ${first_name} was INCORRECT. Corrected in Col 1./" $root\-all-raw3.csv
                            fi
                        fi
                    else
                        sed -i '' -e "${i}s/\$/, Coord search did not find IRAS ID./" $root\-all-raw3.csv
                    fi
                fi
                break
            fi
        done
    #Section 3c finished.
    
    #Section 3d: On the first try (starting with 2 arcseconds), we may get a source that has an IRAS name. If that's the case, we look at the following section. From there -- and again because of formatting, + and - declinations take separate paths -- , we get the alternate names and the Simbad-provided RA and dec (both decimal and sexag). Sometimes, the RA and dec are not given in decimal, in that case we need to take the RA and dec in sexagesimal and go to a converter (I chose http://www.astrouw.edu.pl/~jskowron/ra-dec/ as my converter). Then, sometimes the code might behave weirdly in that the final RA and dec may be blank or contain a value that is not expected, in which case, if they do, the input source is looped again to get a different, more correct value. This serves as a check. Assuming the value is correct, we then replace the input name with the Simbad-provided name, same with RA and dec. We place the RA and dec in sexag and the original in other columns, so we can compare. We also note what the offset between the source from Simbad and the input source and also the search radius (in this case, it's 2 arcseconds). Then, we compare the original and Simbad-provided names. If the Simbad name is IRAS and the original is not, we say that an IRAS name is found. If the Simbad name is an IRAS name, and it matches the original, we say No issues. If the Simbad name is an IRAS name, and the original is also an IRAS name, but the two are different, we say that the original is incorrect, and that we have replaced it.
    elif [[ $name == *IRAS* ]]
    then
        if [[ $dec == \+* ]]
        then
            new_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep  "IRAS [A-Z0-9]" | sed 's/.*"//;s/<.*//;s/>//' | grep  "IRAS [A-Z0-9]" | tail -1)
            new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
            dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)
            altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
            altname="${altname//,/:}"
            
            
            if [[ $new_name != IRAS\ * ]]
            then
                new_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//')
            fi
            
            if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
            then
                dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
            fi

        else
            new_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep  "IRAS [A-Z0-9]" | sed 's/.*"//;s/<.*//;s/>//' | grep  "IRAS [A-Z0-9]" | tail -1)
            new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
            dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1 )
            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1 )
            altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
            altname="${altname//,/:}"
            
            if [[ $new_name != IRAS\ * ]]
            then
                new_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//')
            fi
            
            if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
            then
                dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1 )
                new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1 )
            fi
        fi
        
        if [[ $new_ra == "271.087458" ]] || [[ $dec3 == "-29.519319" ]]
        then
            i=$(( i-1 ))
        elif [[ $new_name == "" ]] || [[ $new_ra == "" ]] || [[ $dec3 == "" ]] || [[ $new_ra_s == "" ]] || [[ $new_dec_s == "" ]] || [[ $altname == "" ]]
        then
            i=$(( i-1 ))
        else
            offset=$(python offset.py $ra $dec $new_ra $dec3)
            offset="${offset//arcsec/}"
            
            sed -i '' -e "${i}s/$first_name/$new_name/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$ra/$new_ra/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$dec/$dec3/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$ra_s/$new_ra_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$dec_s/$new_dec_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $altname/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $ra_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $dec_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $offset/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, 2/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $first_name/g" $root\-all-raw3.csv
            echo "Has IRAS name: replace ${first_name} with ${new_name}, ${ra} with ${new_ra}, ${dec} with ${dec3}"

            if [[ $new_name == *IRAS* ]]
            then
                if [[ $first_name !=  *IRAS* ]]
                then
                    sed -i '' -e "${i}s/\$/, ${first_name} replaced with IRAS in Col 1./" $root\-all-raw3.csv
                else
                    if [[ $first_name == $new_name ]]
                    then
                        sed -i '' -e "${i}s/\$/, No issues./" $root\-all-raw3.csv
                    else
                        sed -i '' -e "${i}s/\$/, ${first_name} was INCORRECT. Corrected in Col 1./" $root\-all-raw3.csv
                    fi
                fi
            else
                sed -i '' -e "${i}s/\$/, Coord search did not find IRAS ID./" $root\-all-raw3.csv
            fi
        fi
    #Section 3d finished.
    
    #Section 3e: On the first try (starting with 2 arcseconds), we may get a source, but it DOESN'T have an IRAS name. If that's the case, we look at the following section. From there -- and again because of formatting, + and - declinations take separate paths -- , we get the alternate names and the Simbad-provided RA and dec (both decimal and sexag). Sometimes, the RA and dec are not given in decimal, in that case we need to take the RA and dec in sexagesimal and go to a converter (I chose http://www.astrouw.edu.pl/~jskowron/ra-dec/ as my converter). Then, sometimes the code might behave weirdly in that the final RA and dec may be blank or contain a value that is not expected, in which case, if they do, the input source is looped again to get a different, more correct value. This serves as a check. Assuming the value is correct, we then replace the input name with the Simbad-provided name, same with RA and dec. We place the RA and dec in sexag and the original in other columns, so we can compare. We also note what the offset between the source from Simbad and the input source and also the search radius (in this case, it's 2 arcseconds). Then, we compare the original and Simbad-provided names. If the Simbad name is not an IRAS name, we say that no IRAS name has been found.
    else
        if [[ $dec == \+* ]]
        then
            new_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//')
            new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/%2b.*//')
            dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*%2b//' | sed 's/./+&/')
            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/+.*//' | tail -1)
            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*+//;s/./+&/' | tail -1)
            altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
            altname="${altname//,/:}"
            
            if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
            then
                dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/' | sed 's/+/%2B/')
                new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/+.*//' | tail -1)
                new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}${new_dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*+//;s/./+&/' | tail -1)
                
            fi
        else
            new_name=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'The astronomical object called' | sed 's/.*called //;s/is.*//')
            new_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/-.*//')
            dec3=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'radius: ' | sed 's/.*c=//;s/&-c.*//' | tail -1 | sed 's/.*-//' | sed 's/./-&/')
            new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/-.*//' | tail -1)
            new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep 'NAME="Coord"' | grep 'ID="Coord"' | sed 's/.*VALUE="//;s/".*//' | sed 's/.*-//;s/./-&/' | tail -1)
            altname=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep '^<A HREF=' | grep 'Dic-Simbad' | sed 's/<\/A>//;s/.*"//;s/<.*//;s/>//' | tr '\n' ';')
            altname="${altname//,/:}"
            
            if [[ $new_ra != [0-9][0-9][0-9][.]* ]]
            then
                dec_ra=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/ //g' | sed 's/./J&/')
                new_ra=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/ .*//')
                dec3=$(curl --silent "https://www.astrouw.edu.pl/~jskowron/ra-dec/?q=${dec_ra}" | grep "d [0-9]" | grep 'href="http://www.astrouw.edu.pl/~jskowron/ra-dec/?q=' | sed 's/.*d //;s/".*//' | tail -2 | head -1 | sed 's/.* //')
                new_ra_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/-.*//' | tail -1)
                new_dec_s=$(curl --silent "https://simbad.u-strasbg.fr/simbad/sim-coo?Coord=${ra}+${dec}&CooFrame=ICRS&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=ICRS-J2000&Radius=2&Radius.unit=arcsec&submit=submit+query&CoordList=" | grep "wavelength class for the origin of the measure" | sed 's/<font color="grey"><i>//g;s/<\/i><\/font>//g;s/<SPAN title="wavelength class for the origin of the measure">//g' | sed 's/.*-//;s/./-&/' | tail -1)
            fi
        fi
        
        if [[ $new_ra == "271.087458" ]] || [[ $dec3 == "-29.519319" ]]
        then
            i=$(( i-1 ))
        elif [[ $new_name == "" ]] || [[ $new_ra == "" ]] || [[ $dec3 == "" ]] || [[ $new_ra_s == "" ]] || [[ $new_dec_s == "" ]] || [[ $altname == "" ]]
        then
            i=$(( i-1 ))
        else
            offset=$(python offset.py $ra $dec $new_ra $dec3)
            offset="${offset//arcsec/}"
            
            sed -i '' -e "${i}s/$first_name/$new_name/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$ra/$new_ra/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$dec/$dec3/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$ra_s/$new_ra_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/$dec_s/$new_dec_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $altname/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $ra_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $dec_s/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $offset/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, 2/g" $root\-all-raw3.csv
            sed -i '' -e "${i}s/\$/, $first_name/g" $root\-all-raw3.csv
            echo "No IRAS name: replace ${first_name} with ${new_name}, ${ra} with ${new_ra}, ${dec} with ${dec3}"

            if [[ $new_name == *IRAS* ]]
            then
                if [[ $first_name !=  *IRAS* ]]
                then
                    sed -i '' -e "${i}s/\$/, ${first_name} replaced with IRAS in Col 1./" $root\-all-raw3.csv
                else
                    if [[ $first_name == $new_name ]]
                    then
                        sed -i '' -e "${i}s/\$/, No issues./" $root\-all-raw3.csv
                    else
                        sed -i '' -e "${i}s/\$/, ${first_name} was INCORRECT. Corrected in Col 1./" $root\-all-raw3.csv
                    fi
                fi
            else
                sed -i '' -e "${i}s/\$/, Coord search did not find IRAS ID./" $root\-all-raw3.csv
            fi
        fi
    fi
        
    i=$(( i+1 ))
    n=$(wc -l < $root\-all-raw3.csv)
done
#Section 3e finished.
#Section 3 finished

#Section 4: search.sh is now done querying Simbad. We sort root-all-raw3.csv by ascending RA again just to make sure. Then, going from root-all-raw3.csv to root-all-final.csv, we sort root-all-raw3.csv by offset (this does not change root-all-raw3.csv itself) then delete duplicates that share RA and dec (deg) because for the duplicates, we want to keep the one that was found at the lowest search radius. Then we sort again by ascending RA and output to root-all-final.csv. For the ones deleted, going from raw3 to final, we output that to root-all-final-rejects.csv. We add the header to the rejects file as well. Then we create the file that picks the sources that were found at less than a user-specified search radius (root-all-final_leThres.csv). We also create another rejects file that contains the stars that WERE NOT picked (root-all-final_leThres_rejects.csv). Then, take the sexagesimal coordinates only from the _leThres file and create the coords file, which will be used when inputting on the VizieR website (root-all-final_leThres_coords.csv). We can then use this list to input into VizieR!

sort -nk2 -t"," -o $root\-all-raw3.csv $root\-all-raw3.csv

sort -nk11 -t"," $root\-all-raw3.csv | awk 'BEGIN{FS=", ";OFS=", "} !seen[$2,$3]++' | sort -nk2 -t"," > $root\-all-final.csv

comm -23 <(sort $root\-all-raw3.csv) <(sort $root\-all-final.csv) > $root\-all-final_rejects.csv

sort -nk2 -t"," -o $root\-all-final_rejects.csv $root\-all-final_rejects.csv

sed -i '' -e '1 i\
#Name, New RA (deg), New Dec (deg), New RA (sexag), New Dec (sexag), Vizier Cat, Ref Number, Alt Names, Orig RA (sexag), Orig Dec (sexag), Offset (arcsec), Search Radius (arcsec), Original Name, Comment on Name' $root\-all-final_rejects.csv

awk -v thresh=$thresh 'BEGIN{FS=", ";OFS=", "} {if($12<=thresh)print $0}' $root\-all-final.csv > $root\-all-final_leThres.csv

comm -23 <(sort $root\-all-final.csv) <(sort $root\-all-final_leThres.csv) | sort -t"," -nk2  > $root\-all-final_leThres_rejects.csv

sed -i '' -e "s/\$/, $thresh/g" $root\-all-final_leThres.csv

sed -i '' -e '1 i\
#Name, New RA (deg), New Dec (deg), New RA (sexag), New Dec (sexag), Vizier Cat, Ref Number, Alt Names, Orig RA (sexag), Orig Dec (sexag), Offset (arcsec), Search Radius (arcsec), Original Name, Comment on Name, Threshold Offset (arcsec)' $root\-all-final_leThres.csv

sed -i '' -e '1d' $root\-all-final_leThres_rejects.csv

sed -i '' -e "s/\$/, $thresh/g" $root\-all-final_leThres_rejects.csv

sed -i '' -e '1 i\
#Name, New RA (deg), New Dec (deg), New RA (sexag), New Dec (sexag), Vizier Cat, Ref Number, Alt Names, Orig RA (sexag), Orig Dec (sexag), Offset (arcsec), Search Radius (arcsec), Original Name, Comment on Name, Threshold Offset (arcsec)' $root\-all-final_leThres_rejects.csv

awk 'BEGIN{FS=", ";OFS=" "}{print $4,$5}' $root\-all-final_leThres.csv | sed '1d' > $root\-all-final_leThres_coords.csv

sed -i '' -e 's/^ *//;s/ *$//;s/  */ /;' $root\-all-final_leThres_coords.csv
#section 4 finished.
