################################
##                            ##
##       Miranda Chang        ##
##                            ##
##  Mentor: Raghvendra Sahai  ##
##                            ##
##       JPL MSP 2021         ##
##                            ##
##                            ##
##                            ##
##                            ##
##   Master Description and   ##
##      Record of Work        ##
##                            ##
################################

======== Contents ========
  0 - Summary
    0.1 - Purpose
    0.2 - Procedures
    0.3 - Directories and Notes
    0.4 - Project Steps and Results
  1 - Unix and Python Notes
    1.1 - Unix
    1.2 - Python
  2 - Programs and Pipeline Scripts
    2.1 - Compiling XX Stars Program
    2.2 - Calculating Offsets Between Input and Output Coordinates Program
    2.3 - Making Gaia EDR3 and Bailer-Jones Data Readable Program
    2.4 - Getting Photometry and Dictionary Program
        2.4.1 - Input file
        2.4.2 - Coordinates file
    2.5 - SED Plotting Program
    2.6 - Color-Color Plotting Program
        2.6.1 - Input Colors file
    2.7 - Querying Images from IRSA
        2.7.1 - Image Parameters file
        2.7.2 - Parameter Records
    2.8 - Making Postage Stamp Mosaics
        2.8.1 - Mosaic Parameters file
    2.9 - Post-RGB Star Candidate List
    2.10 - Creating LaTeX Tables Containing SED Input Data
    2.11 - Creating LaTeX Tables Containing Bailer-Jones/GAIA Data
        2.11.1 - Input Stars Bailer-Jones/GAIA File
  3 - Modeling SEDs of post-RGB stars
  4 - Ongoing/Resolved Issues (TO-DO)
  5 - Miscellaneous Info

======== 0 - Summary ========

      ------- 0.1: Purpose -------
      We continue to discover new types of stars in interstellar space over time. One type that we believe to be a new category is post-red giant branch (RGB) stars. These types of stars are similar to post-asymptotic giant branch (AGB) stars, which are stars that have just evolved off the AGB after having lost most of their stellar envelope via a massive wind and are increasing in temperature while having constant luminosity. What separates post-RGB stars from these types is that post-RGB stars have significantly lower luminosity (100 - 2,500 L☉, much less than 5,000 - 10,000 L☉ from post-AGB stars) and have evolved off the RGB through interacting with another star (Kamath et al. 2016). So far, only one stellar object — the Boomerang Nebula —  can be classified as a post-RGB star because of its low luminosity (300 L☉), which was determined by observing how much it was losing, thus putting an upper limit on how faraway star is and on its luminosity, and how it appears to have resulted from a binary system while the star was still on the RGB (Sahai et al. 2017). Stars in the Large and Small Magellanic Clouds were also found to be strong candidates because of their luminosities that can be verified by the distances to the LMC and SMC (Kamath et al. 2016). 
      
      There are few stars that fall under the classification of post-RGB because it is difficult to find distances to red giant stars. We can reliably measure distances up to 500 parsecs (pc) thanks to the Hipparcos astrometric mission, but most red giant stars are relatively rare, so they lie at distances greater than 500 pc. In addition, AGB stars are very big (approximately 1-2 astronomical units (AU) in diameter), and the brightness of their surface is not uniform and changes with time, which adds noise to the parallax. Furthermore, dust surrounds these stars, making the star faint, and thus, harder to measure the distance to the star. The lack of post-RGB stars form the basis of this project, and the goal is to discover more potential post-RGB stars. While we could study stars in the LMC and SMC, we do not have enough angular resolution to easily study them. It is more beneficial when they are closer to us where we can determine their physical properties, necessary for constraining theoretical models.
      
      Based on what we know about post-RGB stars, we plan to start searching for more by looking at OH/IR stars. These are another type of star that have ejected much of their stellar envelope in a dusty, molecular wind, making these stars bright sources of infrared and OH maser emission. These features occur when a star is in the later stages of AGB evolution. Thus, post-AGB stars make a large portion of these stars. We also believe that post-RGB stars can also be found within the list of OH/IR stars because of their similarities to post-AGB stars. From there, these were the following steps we took (also see procedures for more detailed steps):
      
      - Compile OH/IR stars from VizieR, get more accurate coordinates from Simbad, get stars that have an offset of <= 10 arcseconds (step 1)
      - Query Gaia EDR3 and Bailer-Jones to get trignometric parallaxes and distances for these stars and filtered out stars that did not have them (step 2)
      - Used get-dict-phot.py to get photometry data on these stars (step 3)
      - Used phot_plotter.py to calculate bolometric fluxes and luminosities (step 4)
      - Picked out stars that satisfied the luminosity criteria (> ~30 L☉ and < 1500 L☉) --> these are our post-RGB stars (step 5)
      - Modeled spectral energy diagrams (SEDs) and deteremined mass-loss rates (steps 6 & 7)
      - Created image mosaics from IRSA (step 8)

      ------- 0.2: Procedures -------
      Step 1 - Getting OH/IR stars: 
      Type "OH/IR stars" into VizieR (https://vizier.u-strasbg.fr/viz-bin/VizieR):
          -When viewing the catalogs, some catalogs will be missing information. Also, some catalogs will contain subparts called tables. In order for the catalog to be used, or part (tables) of the catalog, it NEEDS to contain name, right ascension, and declination. To check if the catalog satisfies that condition, click the hyperlink (it starts with J/..). On that page, there will be a submit button, right next to the reset all button. Click the submit button. You should see the name, right ascension, and declination within the first several columns. If one of those is missing for one table, but is visible in other tables, the catalog is good. If all tables/catalog are missing at least one of those features, you may have to do more digging and view the "Read Me" file on the query page to go further and find the original catalog. Otherwise, DO NOT use that catalog. The following catalogs apply (aka DO NOT USE):
          -Amiri+, 2012 (J/A+A/538/A136)
          -Lewis 1997 (J/ApJS/109/489). However, when looking at the Read Me, the stars are actually taken from the Arecibo 1612 mission, so enter "Arecibo 1612" in the search bar on the main page and query those (Lewis 1994 J/ApJS/93/549).
          -Lepine+, 1995 (J/A+A/299/453)
          Also, don't use the following catalogs because they are included in others:
          -Chengalur+ 1993 (J/ApJS/89/189) (included in Engels+ 2015)
          -Chen+, 2001 (J/A+A/368/1006) (included in Engels+ 2015)
          Next, when it comes to right ascension and declination, we need to have the coordinates in both degrees and sexagesimal. With each catalog, and even each table, the coordinates are given in one of those two. The other coordinates need to be supplied by VizieR. So, with each catalog, check if the coordinates are given in sexagesimal (h:m:s) or degrees. After checking, if the coordinates by the catalog are given in sexagesimal, then we need to select decimal for VizieR to supply (for more on this, we will see in the next following steps). Sometimes, the tables within each catalog have coordinates different from each other (meaning some tables will be given in sexagesimal and some in degrees, all within in the same catalog). In that case, you will have to go table by table and download the information. Now, we move on to downloading the information. Download everything (ex: photometry, coordinates) about star, so check all columns before pressing submit. Change 50 to unlimited and HTML table to .tsv file (tab-separated) on the left hand side. For where it says "Position in:", choose the coordinate type that IS NOT given by the catalog. Repeat for remaining catalogs. Rename and place all of these .tsv files into a folder (you can name it whatever you want, but the name of the folder will be the root name for all of the outputs from running search.sh) (location ~/mchang-JPL2021/root). 
          To make things clearer, I will give an example of the above. Let me take the Sevenster catalog. When I see the information displayed about Sevenster, I see that for the first table (J/AJ/123/2772/AGB), the coordinates given are in degrees. For the other two tables, (J/AJ/123/2772/doublede and J/AJ/123/2772/tablevla), the coordinates are given in sexagesimal. So, because of this, I will have to split up the catalog into tables that have degrees and tables that have sexagesimal. So, J/AJ/123/2772/AGB will be one file and J/AJ/123/2772/doublede and J/AJ/123/2772/tablevla will be in another file. For J/AJ/123/2772/AGB, when I'm downloading this file, I will select "sexagesimal" (where it says "Position in:"). Likewise J/AJ/123/2772/doublede and J/AJ/123/2772/tablevla, I will select "decimal."
          So far, the catalogs I have are:
          -J/ApJS/93/549
          -J/AJ/123/2772
          -J/A+AS/92/43
          -J/A+AS/127/185
          -J/A+A/391/967
          -J/A+A/446/773
          -J/AJ/127/501
          -J/A+AS/128/35
          -J/A+A/431/779
          -J/A+A/562/A68
              -by extension of above, J/A+A/562/A68/maps
          -J/MNRAS/505/6051
          There may be more catalogs coming out in the future, so check back to see if any new catalogs pop up. In addition, Dr. Sahai might send some catalogs for you to check in VizieR. Then, run the search.sh file (/mchang-JPL2021/search.sh) (see section 2.1 for more details and for how to run it). In short, to run the program, if it's the first time or if you get an error when you just type ./search.sh root threshold, type chmod +x search.sh, press enter, then ./search.sh root threshold (root and threshold are user-specified). Every time after that, you can just type ./search.sh. Head's up, this program TAKES A VERY LONG TIME TO RUN (>12 hours). The final output from search.sh is root-all-final_leThres_coords.csv is used in Step 2. (If you have a Mac and to make sure your computer is still awake when running the program, type caffeinate ./search.sh root threshold.)
          
------------------------------------------------------------------------------------------------          
          
      Step 2 - Getting Trignometric Parallaxes and Distances: 
      Go to the VizieR catalog website. Type "Gaia eDR3" in the search bar and check the boxes for I/350 and I/352 and press "Show table details." Then, check the boxes "I/350/gaiaedr3"  and "I/352/gedr3dis." Then, select "Query selected Tables." Change the search radius from 2 arcmin to 2 arcsec (or some other specified search radius). Check all columns (the main table). On the left hand side, change 50 to unlimited and HTML table to .tsv file (; separated). For the section "Position in," you want to choose the same type as your input coordinates (ex: if you used sexagesimal, select "sexagesimal"). At the top, click "list of targets" where you upload root-all-final_leThres_coords.csv. After all of that, you can press "submit." When you are finished downloading, give it a name (I chose root-all-final_Thres_eDR3-bailjons.tsv, but again up to you) and put it in the main directory. Then, run edr3-bailerjones-dist.sh (you can see more details in 2.3). The outputs root_coords-master.csv and rootx_coords.csv should be moved to the same directory as get-dict-phot.py, phot_plotter.py, and color-color_plotter.py (you can do this by using the mv command or if it's really difficult, dragging and dropping the files). rootx_coords.csv will be used in step 3.
      
------------------------------------------------------------------------------------------------ 
      
      Step 3 - Get Photometry data from VizieR: 
      In step 2 (edr3-bailerjones-dist.sh), I had split up the coords-master file into x number of files. This is because if I were to run a humongous list of stars, there might some issues with querying VizieR in that sometimes querying a star in VizieR might cause it to return a ValueError. Therefore, I broke it up into x number of files (in my case I broke it up into 20 different files containing roughly 75 lines each), so in case VizieR returns an error, it would be easier to find the error in a smaller list of stars compared to one giant list. If something magical happens where VizieR returns no errors with the giant list, then there is no need to break the file up, and we can rename the root_coords-master.csv file as root_coords.csv and run get-dict-phot.py. However, in the meantime, we will stick with our broken up list. In this case, make x number of copies of the .inp file, make sure they have the same root name as their coords file. Run get-dict-phot.py for each coords file (see section 2.4). HEADS UP: the coords.csv contains coordinates from Gaia EDR3 (2020). Choose a search radius in arcsec and specify it in the .inp file. See section 2.4.1 for the .inp file for querying the necessary catalogs. Enter in terminal: python get-dict-phot.py root.inp. 
      Now what if you run into ValueErrors? As soon as you get a ValueError, see what star that is in the termianl window. Then go to the corresponding coords.csv and find that star. Move that star and its information to the top of the coords file and run get-dict-phot.py again (this will make the affected star the first to run). If a ValueError is still returned (you will found out pretty quickly), then comment out that line in the coords.csv by putting a # in front of the line. Run get-dict-phot.py on this coords file again. Repeat this for the rest of ValueErrors you get. At the end of running all the coords file, grab all the lines that start with # from each coords file (you can do grep '#') and place them into another coords file caused root-leftovers_coords.csv. Run get-dict-phot.py on this file. It is possible that they may start working again as ValueErrors tend to change every few hours or so (in my case, most of my leftover stars worked again!). The files in the outputs from get-dict-phot.py - folder rootx-yasec - will be used in step 4 (i.e. _phot.csv, _spec.csv, _integresults.txt, etc.). The sed files (really just the .dat) in sed-for-MoD in rootx-yasec will be used in step 7.
      
------------------------------------------------------------------------------------------------ 
      
      Step 4 - Calculate Bolometric Fluxes and Luminosities using Distances from Step 2: 
      Run phot_plotter.py for each coords file (see section 2.5). This outputs the bolometric fluxes and calculated luminosities, along with info from the Gaia EDR3 and Bailer-Jones .tsv. Make sure get-dict-phot.py is run first as we need the outputs from that program for this program. Enter in terminal: python phot_plotter.py root.inp. The dist-lum.csv in each rootx-yasec folder are used in step 5.
      
------------------------------------------------------------------------------------------------ 
      
      Step 5 - Find post-RGB Star Candidates and Criteria:
      To get the post-RGB candidate list, we select stars that have an actual luminosity of > ~30 L☉ and < 1500 L☉. Run the program prgb_candidates.sh (type chmod +x prgb_candidates.sh, hit enter, then ./prgb_candidates.sh root, see section 2.9 for more info). Note that you may have to alter the code yourself depending on how you broke up your master coords (i.e. how many files you have, did you give different names to the files, etc.). This outputs stars that match the above criteria and also appends their OH maser line velocities.
    
------------------------------------------------------------------------------------------------ 

      Step 6 - Creating Image Mosaics for Stars
      Then, to get postage stamp images for our candidates, run the program postage_stamp.py using the following command: python postage_stamp.py ohir-img.inp (see section 2.7 for further information). It might be helpful to have ds9 installed as it will allow you to look at the fit files and troubleshoot any problems pertaining to them. To install, type the following in the terminal: brew install --cask saoimageds9. Once you are finished running postage_stamp.py, run the next program,  which is mosaic.py (this creates the mosaic layout and includes the outputs from postage_stamp.py). To run, type: python mosaic.py ohir-mosaic.inp (see section 2.8 for further information). If it is your first time running mosaic.py, set the usr_flg to 0 in the ohir-mosaic.inp file so it will create the .csv's containing the default min and max intensities, to which you will edit later and thus edit the noise in the mosaic images. 
      In this step, you may need to pip install regions and pip install cvskit if you get any errors related to these two.
 
------------------------------------------------------------------------------------------------ 

      Step 7 - Setting up MoD: 
      We are going to make models of the post-RGB candidates. Follow the series of steps to download the MoD program:
      1) go to http://homepage.oma.be/marting/codes.html
      2) download the two versions of this code, which you will need to make models of the pRGB
        candidates that we find:
            
            "The latest (July 2018) External release 6 (Internal release 16a) is available in this tar file"
            "A test version: Release 7 (November 2020) is available in this tar file"

         when you untar, you will get a dir called MoD...

      3) install* the test version using the source code MyDusty_v16a. (as per directions in
        MoDManual) in the subdir oh26.5 -- you will run into some problems because there are
        references to directory paths that are not valid for you, as these contain:

        /home/marting

        you simply need to globally change /home/marting with the full pathname to the MoD
        directory. e.g., on my computer, I have /home/sahai/prog/MoD-V7 as the main program
        directory, so I replaced /home/marting with /home/sahai/prog/MoD-V7.
        To do this, you will use an example of this command: 
        
sed -i '' -e 's/home\/marting/\/Users\/miranda\/mchang-JPL2021\/V7/g' ~/mchang-JPL2021/V7/oh26.5/MyDusty_v16a.f 
        
        NOTE that the above command is for Mac OS (It requires -e). Most likely other operating systems will NOT use '' -e. Also, Fortran has a character limit for each line so make sure your path name is not too long. You will get an error if it is when you run the gfortran command seen below. If it is too long, you might want to move the MoD to another directory, so the path name is not as long.

        the code requires that you have the plotting package pgplot installed as well. this is
        available at: https://sites.astro.caltech.edu/~tjp/pgplot/

        if you run into problems that you can't resolve, you can also install MyDusty_v16a_noPGP.f, which is a version of the MyDusty code that does not call the pgplot package.
        
        Trying to install PGPlot didn't work for me. I ran into a roadblock with one of the steps. I used the following steps:
        
        sudo mkdir /usr/local/pgplot
        sudo mkdir /usr/local/src/pgplot
        sudo chown miranda /usr/local/pgplot
        sudo chown miranda /usr/local/src/pgplot
        tar -xf ~/Downloads/pgplot5.2.tar.gz -C /usr/local/src
        cd /usr/local/pgplot
        cp /usr/local/src/pgplot/drivers.list /usr/local/pgplot
        vi drivers.list
        Go over to the ! for all the GIDRIV and PSDRIV lines. Make sure cursor is right on the exclamation point and press "x" to delete the exclamation points (that deletes them and uncomments the line). We need those lines for plotting.
        After deleting the necessary !, Type in ":". You will see the : appear in the bottom left of the window. Then type in "wq" and then hit enter. That saves the changes you made to the doc.

        /usr/local/src/pgplot/makemake   /usr/local/src/pgplot  Linux g77_gcc (last step I tried)
        
        I couldn't get past the last step. I tried to install via XCode too, but I didn't have enough room on my computer to install the newest operating system. Therefore, I had no choice but to install no PGPlot, but that is ok as there is a way to model without PGPlot. If you still need help, contact William Adler for help.
        
        *To install the test version, if you have the pgplot installed, use the command: 
        gfortran  -fbounds-check -O -o MyDusty.exe  MyDusty_v16a.f -lpgplot -std=legacy
        If you DON'T have pgplot installed, use the command:
        gfortran  -fbounds-check -O -o MyDusty.exe  MyDusty_v16a_noPGP.f -std=legacy
        
        You will have to input the file path for MyDusty.exe and MyDusty_v16a.f/MyDusty_v16a_noPGP.f to match yours (if needed) (ex: for me, MyDusty.exe has the path ~/mchang-JPL2021/V7/oh26.5/MyDusty.exe)
        
        To check if the version you installed is working, type one of the following into terminal: ./MyDusty.exe oh26.5_in. If it runs through with no issues or errors popping up, then it works (it will take maybe 45 minutes for the program to run).
        
        To plot data, enter this command: gnuplot -p PlotMyDusty_oh003.234.gnu (replace oh003.234 with the root) (if you don't have gnuplot, type brew install gnuplot)
        
------------------------------------------------------------------------------------------------ 
        
      Step 8 - Setting Up a Directory for Five of Your Post-RGB Stars in MoD: 
      With MoD set up and the output from get-dict-phot.py ready (specifically the sed2s), we can start modeling. Before we start, choose the five best candidates using the following criteria:
      
      1) postage stamp image mosaic unambiguously shows that candidate star is the same as the one ID'd in Bailer-Jones/eDR3. so the circle should be centered on the bright star in the images, and there should not be very close sources nearby (within an offset of ~1 arcsec from  center of circle) which could cause confusion
      2) most extensive photometry (so MoD model results better constrained)
      3) OH/IR velocities available
      
      Once you are done selecting, copy all the files of oh26.5 into another directory that has the name of the star that you are testing with (you can do cp -r root1 root2). To change the root name of the files in the folder, do brew install rename and then rename 's/root1/root2/' *. 
      Dr. Sahai gave me clearoutputfiles.sh, which clears all the outputs of that folder, so run it. Make sure you have your sed2_name.dat file in your new directory (it is from step 3 in the rootx-yasec folder in the folder sed-for-MoD). Use the following code to copy it over: cp ~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/sed-for-MoD/sed2_name.dat sed2_root.dat . Then, within that directory, go to the file that ends with _in. There you will see a list of parameters (it also might be handy to have the MoD manual with you - when you download MoD, it might be there as a .pdf, but for me, I had to use Overleaf and compile MoDManual.tex) that will plot the best fit for the data. Every star is different, so play around the values and pick the ones that best fit (see section 5 for more details about parameters). 
      After selecting the parameters, run the commands mentioned in the previous step (./MyDusty.exe root_in -> gnuplot -p PlotMyDusty_root.gnu.) For the plot, you need to have rted_root.dat and sed2_root_nufnu.dat. rted_root.dat is created from running ./MyDusty.exe while sed2_root_nufnu.dat you must do manually. It takes the wavelength, nuFnu, and nuFnu error from the sed files. Run the following code to produce the sed2_root_nufnu.dat --> awk -F, '$6 != "  --"' ~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/sed-for-MoD/sed2_name_raw5.csv | awk 'BEGIN{FS=", ";OFS=" "} {print $8,$6,$7}' | sed '1,2d' > sed2_root_nufnu.dat . In PlotMyDusty_root.gnu, make sure you have the correct root names within the file as that changes with each star you have! 
      This might take a while to get right as you view the plots. Make sure you record each of trials and what went right/wrong in your master.inf or another .inf file as it may get long! With every new try, do the following steps:
      cp all the contents of the first trial root folder into  the second trial root folder:
            cp -r root_1  root_2

      then do source clearoutputfiles.sh
        (this removes ALL the OUTPUT files in root_1)

      edit the _in file to change parameters

      run MoD on the new _in file, after its done

      run the *gnu file to make plot
      
      If you find yourself stuck on getting the best fit, you may have to choose a new star.
      
------------------------------------------------------------------------------------------------ 
      
      Step 9 - Creating Figures for Written Reports
      
      After you have gotten the SEDs and mosaics for your select stars, we now create LaTeX tables to help make information about the stars readable should we want to write a report about our findings. Specifically, we take the sed2_name_raw5.csv that was created in step 3 and output the information in raw5.csv into a LaTeX table using the Python program pht_mod.py. See section 2.10 for more information on how to run it. We will also use information from the Bailer-Jones/GAIA file  (ohir-all-final_Thres_eDR3-bailjons-cln.csv) created in Step 2 to create another LaTeX table using the Python program bailjons_motion_data.py. See section 2.11 for more information on how to run it.
      
------------------------------------------------------------------------------------------------ 

      ------- 0.3: Directories and Notes -------
      The terminal instructions for running each program is written in comments at the top of the program.
      
      mchang-JPL2021 (main project dir)
          -- 2mass (creating color-color plots based on 2mass data)
              -- python (python programs from querying from 2mass catalog from a past student)
          -- ohir (raw .tsv files from VizieR)
              -- ohir-clean (the name, ra, dec, catalog name, and ref num from each VizieR catalog)
              -- ohir-clean-v2 (same as ohir-clean but associated with v2 of search.sh being run, see v2 files in past-ohir)
              -- ohir-clean-v3 (same as ohir-clean-v2 but with v3 instead)
          -- python_phot-lum-colors
              -- photometry (python programs to extract photometry data from VizieR)
                  -- ohir-samp_colors (output from color-color_plotter.py)
                  (the following are outputs from get-dict-phot.py and phot_plotter.py)
                  -- ohir01-2.5asec
                  -- ohir02-2.5asec
                  -- ohir03-2.5asec
                  -- ohir04-2.5asec
                  -- ohir05-2.5asec
                  -- ohir06-2.5asec
                  -- ohir07-2.5asec
                  -- ohir08-2.5asec
                  -- ohir09-2.5asec
                  -- ohir10-2.5asec
                  -- ohir11-2.5asec
                  -- ohir12-2.5asec
                  -- ohir13-2.5asec
                  -- ohir14-2.5asec
                  -- ohir15-2.5asec
                  -- ohir16-2.5asec
                  -- ohir17-2.5asec
                  -- ohir18-2.5asec
                  -- ohir19-2.5asec
                  -- ohir20-2.5asec
                  -- ohir-leftovers-2.5asec
                  -- ohir-more-leftovers-2.5asec
                  -- ohir-samp-2.5asec (containing an example of an output from get-dict-phot.py)
                  -- ohir-samp-2.5asec-v0 (containing an example of an output from get-dict-phot.py)
                  -- ohir_img (output from postage_stamp.py)
              -- examples (how to structure your .inps)
          -- past-ohir-all (past versions of ohir-raw outputs)
          -- past-search (past versions of search.sh)
          -- ref-papers (reference papers used for the background of this project)
          
      ------- 0.4: Project Steps and Results -------
      1. Collected lists of OH/IR stars manually from VizieR (~/mchang-JPL2021/ohir)
      
      2. Ran program (~/mchang-JPL2021/search.sh) which:
      -Merged the lists of OH/IR stars into one list (~/mchang-JPL2021/ohir-all-raw1.csv)
      -Sorted the ohir-all-raw1.csv by ascending RA (~/mchang-JPL2021/ohir-all-raw2.csv)
      -Deleted duplicates, obtained alternate names, new coordinates from Simbad and search radius new coordinates were found at, calculated the offset between Simbad and VizieR coordinates (~/mchang-JPL2021/ohir-all-raw3.csv)
      -Deleted duplicates in ohir-all-raw3.csv (~/mchang-JPL2021/ohir-all-final.csv)
      -Sent the deleted duplicates in ohir-all-raw3.csv to another file (~/mchang-JPL2021/ohir-all-final_rejects.csv)
      -Selected stars that had an offset of less than or equal to 10 arcseconds between the Simbad and VizieR coordinates (~/mchang-JPL2021/ohir-all-final_leThres.csv)
      -Sent the stars that had an offset of GREATER than 10 arcseconds into another file (~/mchang-JPL2021/ohir-all-final_leThres_rejects.csv)
      -Extracted just the RA and dec of the stars in ohir-all-final_leThres.csv (~/mchang-JPL2021/ohir-all-final_leThres_coords.csv)
      
      3a. Inputted ohir-all-final_leThres_coords.csv to get Gaia eDR3 and Bailer-Jones data in VizieR (~/mchang-JPL2021/ohir-all-final_Thres_eDR3-bailjons.tsv)
      
      3b. Ran program (~/mchang-JPL2021/edr3-bailerjones-dist.sh), which:
      -Retrieved the headers from the catalogs selected and put them on one line (~/mchang-JPL2021/root-all-final_Thres_eDR3-bailjons-hdr.csv)
      -Outputted ohir-all-final_Thres_eDR3-bailjons.tsv into .csv format along with name and original coordinates of star (~/mchang-JPL2021/ohir-all-final_Thres_eDR3-bailjons-cln.csv)
      -Grabbed stars that have values for the rgeo (distance) column in ohir-all-final_Thres_eDR3-bailjons-cln.csv along with their names, RA and dec in sexagesimal, rgeo in kpc (along with lower and upper limits) proper motion data, and noise data (~/mchang-JPL2021/ohir_coords-master.csv, later moved ~/mchang-JPL2021/python_phot-lum-colors/photometry in next step)
      -Broke up ohir_coords-master.csv into 20 files (~/mchang-JPL2021/ohirx_coords.csv, x means number 1-20, later moved to ~/mchang-JPL2021/python_phot-lum-colors/photometry in next step)
      
      3c. Moved ohir_coords-master.csv and ohirx_coords.csv to photometry folder (~/mchang-JPL2021/python_phot-lum-colors/photometry)
      
      4. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/get-dict-phot.py) using ohirx_coords.csv and ohirx.inp (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx.inp) which retrieved photometry data from VizieR and outputted:
      -SED data, containing fluxes and filters (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/sed-for-MoD)
      -Rejects (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_rejects.txt)
      -What columns that data are in (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_hdr_column.csv)
      -What coordinates were used (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx_coords.csv)
      -What missions were used (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_missions.csv)
      -Photometry data (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_phot.csv)
      - ~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_phot_raw.csv
      
      5. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/phot_plotter.py) using ohirx_coords.csv and ohirx.inp (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx.inp) which calculated bolometric flux and luminosities and outputted:
      -SED plots (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/name_fnu-wav.pdf, ~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/name_nufnu-wav.pf, and ~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/name_spectrum.pdf)
      -.txt file with fluxes (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_integresults.txt)
      -.csv file with velocities, fluxes, and luminosities (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/ohirx-2.5asec_dist_lum.csv)
      
      6. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/prgb_candidates.sh) which took stars in ohirx-2.5asec_dist_lum.csv and outputted those that are > ~30 L☉ and < 1500 L☉ along with OH maser velocities found in ohirx-2.5asec_phot.csv (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir-prgb-final-candidates.csv)
      
      7. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/postage_stamp.py) which got IRSA images for stars in ohir-prgb-final-candidates.csv using ohir-img.inp (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir-img.inp) --> (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir_img/name_imgs)
      
      8. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/mosaic.py) which created mosaics using the images (~/mchang-JPL2021/ohir_img/name_imgs) and outputted luminosity, velocity, and distance info (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir-prgb-final-candidates.csv). Also used ohir-mosaic.inp (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir-mosaic.inp)
      -Mosaic (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir_img/name_imgs/name_ra_dec_mosaic.png)
      -.csv containing filter name and min and max intensities of each filter (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohir_img/name_imgs/ohir_ra_dec-min-max.csv)
      
      9. Ran MoD (~/V7/root_iteration) which created SED models and outputted fitted characteristics, such as temperature, luminosity, of star (~/V7/root_iteration/root.inp) using filter names along with fluxes that were created in get-dict-phot.py (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/sed-for-MoD/sed2_name.dat), guess parameters (~/V7/root_iteration/root_in), and wavelengths, nuFnu, and nuFnu error of filters (~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/sed-for-MoD/sed2_name_raw5.csv)
      
      10. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/pht_mod.py) which created LaTeX table containing filter name, wavelength, flux, flux error, offset, and name of catalog (~/mchang-JPL2021/python_phot-lum-colors/photometry/name_pht_MoD.tex) using ~/mchang-JPL2021/python_phot-lum-colors/photometry/ohirx-2.5asec/sed-for-MoD/sed2_name_raw5.csv.
      
      11. Ran program (~/mchang-JPL2021/python_phot-lum-colors/photometry/bailjons_motion_data.py) which created LaTeX table containing info from Bailer-Jones/GAIA, specifically name of source, input RA, input dec, output ra, output dec, offset, parallax, parallax error, proper motion, RA proper motion, RA proper motion error, dec proper motion, dec proper motion error, distance, lower limit distance, and upper limit distance of select stars. The following files were used to create the table: ~/mchang-JPL2021/python_phot-lum-colors/photometry/bailjons_motion_stars.inp and ~/mchang-JPL2021/ohir-all-final_Thres_eDR3-bailjons-cln.csv. The output file is ~/mchang-JPL2021/python_phot-lum-colors/photometry/bailjons_table.tex

======== 1: Linux and Python Notes ========

      ------- 1.1 - Linux -------
      1: To move files, use mv oldDirectory/fileName newDirectory
      2: To rename files, use mv oldName newName (from within file directory)
      3: To list files without a certain extension, use ls -I '*.ext'
      4: To install Python modules via the command prompt, use python -m pip install [package]
      5: To create a .sh file:
          -Type #!/bin/bash at the top of every .sh file (can be created using a text editor) so contents of file can be read as command lines. Then put whatever command lines below that.
          -If you are running a script for the first time, enter chmod +x filename.sh in the terminal to activate shell script.
          -Once activated, every time after that, enter ./filename.sh in terminal after the above step to run your script.
      6: To tar gz a directory, use tar -czvf newfilename.tar.gz filename (outside of directory)
      7: To get the first x lines in a file:
          Ex: head -5 ohir-all.csv | awk 'BEGIN{FS=","}{print $1}' > ~/mchang-JPL2021/2mass/sample.txt
         This prints out just the names of first 5 stars of the master list ohir-all.csv onto a file.
      8: To replace words within file, example (if you are on a Mac): sed -i '' -e 's/home\/marting/mchang-JPL2021\/V7/g' ~/mchang-JPL2021/V7/oh26.5/MyDusty_v16a.f
      9: -std=legacy changes errors to warnings and thus makes the program run.
      10: To install via homebrew, type: brew install program. To uninstall, brew uninstall program.
      11: If terminal gets stuck on command, use (at least for Mac OS) control C to stop program/command running
      12: To output the HTML of a webpage and search certain words off of it: Ex: curl --silent "https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Be+Cen&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id" | grep "IRAS [0-9]" | sed 's/.*"//;s/<.*//;s/>//'
      13: To delete a directory, use rm -rf directory
          

      ------- 1.2 - Python -------
      We are using python3 for all programs.
      1: To install modules from the command prompt (either locally or on server) try:
          python -m pip install [package] *Add on --user if permission issues occur*
         To uninstall you might be lucky and just do pip uninstall package. However, you might have to take a more complicated route if that gives you an error:
         Go to Terminal and type pip freeze to see what programs you have installed. Find the name of the program that you want to uninstall. Then type: pip freeze | grep 'name of program' > requirements.txt. If you have more you would like to uninstall type the above again but use >> requirements.txt instead of > requirements.txt to add on to the existing file. Then pip uninstall -r requirements.txt. You can then delete requirements.txt.


======== 2 - Programs and Pipeline Scripts ========
    ------ 2.1 Compiling XX Stars Program -----
    * Location: /mchang-JPL2021/search.sh *
    - Program Purpose: This code runs through all the files in the 'root' folder, which contains all the various catalogs of root stars downloaded from VizieR, and extracts the name of the star, RA (both deg and sexag), dec (both deg and sexag), name of the catalog, and code of the catalog from each file. A new file is created for each catalog and then stored in the subdirectory 'root-clean'. Then, we concatenate all the files and query the stars in Simbad to get even more accurate coordinates, alternate names, and the differences between the input and output coordinates. Then, we select stars that are below a certain offset. See below for outputs and what differentiate them.
    -Outputs: ~/mchang-JPL2021/root-all-raw1.csv - master list, not sorted 
               ~/mchang-JPL2021/root-all-raw2.csv - master list sorted by RA
               ~/mchang-JPL2021/root-all-raw3.csv - master list (identical duplicates of name, RA, and dec are deleted from raw2), star coordinates are entered into Simbad, better coordinates (both decimal and sexagesimal) are returned, non-IRAS named stars replaced with IRAS names(if possible), alternate names are sought, the search radius of Simbad, and the offset between input and output coordinates are calculated. What this code does: the coordinates are entered into Simbad starting at a search radius of 2 arcseconds. If an object is found, the name of the object (IRAS name, if not, the main name listed on Simbad in bold), the coordinates (sexagesimal and deg) from Simbad, alternate names, the offset, and search radius are printed. If no object is found, we increase  the search radius by 1 arcsec until an object is found.
              ~/mchang-JPL2021/root-all-final.csv - After root-all-raw3.csv, we will delete stars that are identical duplicates in name, RA (deg), and dec (deg)
               ~/mchang-JPL2021/root-all-final_rejects.csv - The stars that were deleted from $root\-all-raw3.csv will be outputted to this file
               ~/mchang-JPL2021/root-all-final_leTHRES.csv - Extracting stars that were found less than number arcseconds from the given coordinates. cut & paste line
               ~/mchang-JPL2021/root-all-final_leTHRES_rejects.csv - The stars that were not chosen from root-all-final.csv will be outputted to this file
               ~/mchang-JPL2021/root-all-final_leTHRES_coords.csv - Extracting just the columns that have the new RA and dec in sexagesimal coordinates (this can be inputted in VizieR)
#Every output acts on the previous (ex: root-all-raw2.csv acts on root-all-raw1.csv, root-all-raw3.csv acts on root-all-raw2.csv)
    - Running Procedure: type ./search.sh root threshold into terminal (if it's the first time or if just typing ./search.sh root threshold gives a permission denied error, type chmod +x search.sh, hit enter, then ./search.sh root threshold). Threshold is in arcsec. If you are worried about your computer falling asleep while the code is running, then type: caffeinate ./search.sh root threshold
    - Ex of running: I want my root name to be 'abc', and my threshold to be 10 arcseconds. In the terminal, I will type ./search.sh abc 10 (or caffeinate ./search.sh abc 10).
    
   ------- 2.2: Calculating Offsets Between Input and Output Coordinates -------
        * Location: ~/mchang-JPL2021/offset.py *
        - Program Purpose: Calculates the offset in arcseconds between two coordinates given of the same star. Right ascension and declination of stars are needed.
        - Input: Right ascension and declination of stars are needed.
        - Output: Returns the offset in arcseconds
        - Running Procedure: use this command: python offset.py RA1 Dec1 RA2 Dec2. NOTE: This line is in search.sh, so when you run search.sh, this command will run automatically.
        
    ------ 2.3 Making Gaia EDR3 and Bailer-Jones Data Readable Program -----
    * Location: /mchang-JPL2021/edr3-bailerjones-dist.sh *
    - Program Purpose: This takes the .tsv output from Vizier and makes the data a readable .csv. It then takes stars that have distance values (rpgeo) in the Bailer-Jones catalog and outputs their name, ra, dec, dist(kpc), dist_lower(kpc), dist_higher(kpc), PM, pmRA, e_pmRA, pmDE, e_pmDE, astrometric_excess_noise, astrometric_excess_noise_sig, ipd_gof_harmonic_amplitude, ruwe 
    - Outputs: ~/mchang-JPL2021/root-all-final_Thres_eDR3-bailjons-hdr.csv - the headers of the data from EDR3 and Bailer-Jones printed onto one line
               ~/mchang-JPL2021/root-all-final_Thres_eDR3-bailjons-cln.csv - The name, original coordinates, and data from .tsv outputted into a .csv (EDR3 and Bailer-Jones data are on one line). Missing data is replaced with -999. If star has multiple data from a catalog, we pair the lines that have the same offset (ex: a star may have two entries of data from EDR3, one data has an offset of 0.4 while the other has an offset of 1.3. There is one entry from Bailer-Jones which has an offset of 0.4. In this case, we put the data that have 0.4 together on one line while the 1.3 data goes onto another line).
               ~/mchang-JPL2021/root_coords-master.csv - grabs lines that have values for the rgeo (distance) column. Those lines (their names, RA and dec in sexagesimal, rgeo in kpc (along with lower and upper limits), proper motion data, and noise data) will be outputted.
               ~/mchang-JPL2021/rootx_coords.csv - breaks up the coords-master.csv output into a specified number (x) of files because of ValueErrors returned by get-dict-phot.py when running large files. Right now, I have the program set up so it breaks coords-master.csv into files of 75 lines (so for me that's around 20 files). If you would like to change this, find the line num_lines=75 (near tthe end of the file) and change it. The outputted files are named from root01 to rootx. These outputs are used for get-dict-phot.py and phot_plotter.py.
    - Running Procedure: type ./edr3-bailerjones-dist.sh root into terminal (if it's the first time or if just typing ./edr3-bailerjones.sh root gives a permission denied error, type chmod +x edr3-bailerjones-dist.sh, hit enter, then ./edr3-bailerjones-dist.sh root).
      Ex for running: if I want my root to be 'abc,' then I will type ./edr3-bailerjones-dist.sh abc
    
      ------- 2.4: Getting Photometry and Dictionary Program -------
        * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/get-dict-phot.py *
        root: root-#asec
        directory: root-#asec
       - Program Purpose: Fetches photometry data from Vizier and IRSA databases and compiles it into a dictionary readable file. Reads molecular line data and spectra data and adds to file. This file (along with others) will appear in a directory which is created inside the directory where the program is located. This directory will be named the root name. In addition, with the available photometry data, for each star, a .dat file along with some raw .csvs are stored in sed-for-MoD. Each file contains the name of the filter, flux, flux error, wavelength, offset, and catalog the data was retrieved from. This is to be used for MoD.
       - Outputs: sed-for-MoD, root_rejects.txt, root_hdr_column.csv, root_coords.csv, root_missions.csv, root_phot.csv, root_phot_raw.csv, XML files if flag is activated, IRAS LRS spectra data if available, Optical spectra if available
       - Running Procedure: First, make sure you run edr3-bailerjones-dist.sh to get root_coords.csv. You technically only really need the name,RA,dec for this program to work but root_coords.csv contains other info that will be useful for phot_plotter.py
       Next, create a file called root.inp in the same directory as the program. This .inp file is formatted as written in section 3.1.
       Run the program from the terminal as: python get-dict-phot.py root.inp
       
            ------- 2.4.1: Input File -------
            This file is named 'root.inp' and is a text file that specifies all the input parameters to run get-dict-phot.py and phot_plotter.py. The format of this file should be as so:
                #all comments are preceded by a hashtag.
                #keywords will be followed by an '=' and then the value of the parameter.
                #sep = string which is the separator of name, ra, and dec in coords file
                #c_format = d or s for decimal or sexagesimal in coord file
                #radius = radius of the search in arcsec
                #there are flags for saving xml files (xmlflag) and plotting spectra (spectrum flag) that are 1 for yes, 0 for no
                radius = 5
                xmlflag = 0
                spectrumflag = 1
                #factors for fluxes
                fluxdiff_werr = 2
                fluxdiff_frac = 0.05
                #input coord error
                inp_cerr = ##
                #astrometric accuracies of missions, NAME_amac = amac in asec
                ...
                #source = xmin,xmax,ymin,xmax, if wanting to scale SED plot axes
                #file formats to save plots as, FILE FORMAT = y or n
                .pdf = y
                .jpg = n
                .png = n
                .eps = n
                distance = distance in kpc, if one distance to be used for ALL objects
                dist_each_source = 0 or 1 for whether distances are provided in the root_coords.csv
                #freq1_GHz = min freq in dist_lum file
                #freq2_GHz = max freq in dist_lum file
                freq1_GHz = 1 (low radio freq)
                freq2_GHz = 345 (high radio freq)
                #imageflag = 0 or 1 for whether to save images of object
                imageflag = 0 or 1
                #hips_XXX = 0 or 1, XXX is name of filter w image
                hips_XXX = 0 or 1
                #spfg = 1 or 0 if additional column in coords file has object type labels
                #mos_XXX = ##, XXX is type of object, ## is number of objects in mosaic
                #row = #, # is # mosaics in row
                #col = #, # is # mosaics in col
                massivestar = 0 or 1 for whether qout.dat is provided (for massive star lists)
                
            ------- 2.4.2: Coordinates File -------
              This is the file used above in the Input file and the get-dict-phot.py file. It technically only really needs the following (should be formatted as):
                #NAME,RA,DEC
                nam1,ra1,dec1
                ...
              where all lines precended by a # are read as comments.

              However, for phot_plotter.py, there is some extra info that may be needed when outputting dist-lum.csv. The file would be formatted as this (if you run edr3-bailerjones-dist.sh, the following will be outputted):
              nam1,ra1,dec1,dist(kpc),dist_lower(kpc),dist_higher(kpc),PM,pmRA,e_pmRA,pmDE,e_pmDE,astrometric_excess_noise,astrometric_excess_noise_sig,ipd_gof_harmonic_amplitude,ruwe

              The extra information can be found in the Gaia EDR3 and Bailer-Jones .tsv. The coordinates used in this file are taken from Gaia EDR3.
       
    ------- 2.5: SED Plotting Program -------
          * Location: ~/mchang-JPL2021/python_phot-lum-colors/photometry/phot_plotter.py *
            root: root-#asec
            directory: root-#asec
           - Program Purpose: Plots SEDs from the data file created in get-dict-phot.py. Plots spectra on top of SEDs and separately if desired. Calculates bolometric flux, distances (This program can obtain distances three ways to calculate luminosities: 1 - calculate distances using radial velocities (in root.inp dist_each_source = 0, distance = 0), 2 - one distance is provided for all sources (dist_each_source = 0, distance = whatever you choose), 3 - distances are provided to each star in the root_coords.csv. (dist_each_source = 1)), and luminosities.
           - Outputs: SED plots, .txt file with fluxes, .csv file with velocities, fluxes, and luminosities
           - Runnning Procedure: In the root.inp file, add the following lines:
              spectrumflag = 0 or 1 depending on whether to plot the spectra separately
              #file format - select y or n depending on whether or not to save the graphs as these formats. At least one must be 'y' in order for the graphs to be saved.
              .pdf = y or n
              .jpg = y or n
              .png = y or n
              .eps = y or n
              spectrumflag = 1 or 0
              freq1_GHz = 1 (low radio freq)
              freq2_GHz = 345 (high radio freq)
              imageflag = 0 or 1
              hips_XXX = 0 or 1
              dist_each_source = 0 or 1 for whether distances are provided in the root_coords.csv
              massivestar = 0 or 1 for whether qout.dat is provided (for massive star lists)
              ...
            You may need to install 'kd' if you get a 'ModuleNotFoundError' in which case use this: pip install git+https://github.com/tvwenger/kd.git
            Same thing for 'pyqt_fit': pip install git+https://github.com/tvwenger/pyqt-fit.git
            To uninstall kd: do, pip freeze | grep 'kd' > requirements.txt then pip uninstall -r requirements.txt
            To uninstall pyqt_fit: pip freeze | grep 'pyqt-fit' > requirements.txt then pip uninstall -r requirements.txt
            Press y after it asks for confirmation
            Run this program from the same directory as get-dict-phot.py. This will be run from the terminal as: python phot_plotter.py root.inp

            (SAHAI to MIRANDA) Ignore the distances and luminosities computed here -- you will modify the program to use distances from eDR3 and/or Bailer-Jones
       
    ------- 2.6: Color-Color Plotting Program -------
      * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/color-color_plotter.py *
        root: root
        directory: root_colors
       - Program Purpose: Plots 2D and 3D color-color diagrams. Then, runs k-means clustering on these plots and labels centroids. Returns a .txt file with the locations of the centroids.
       -This program can take multiple input photometry files and plot their color-color points on the same diagram. Just add more coords files to the root_colors.inp file.
       - Outputs: color-color plots, k-means plots, centroid locations .txt file
       - Running Procedure: Create a file called root_colors.inp in the same directory as the program. This .inp file is formatted as written below (2.2.1). In addition, you may need to modify the code, specifically where the code is looking through a directory creaetd from get-dict-phot.py You may need to change the name of the directory in the code with the search radius that you specified.
        Run this program from the same directory as get-dict-phot.py. This will be run from the terminal as: python color-color_plotter.py root_colors.inp
        
               ------- 2.6.1: Input Colors File -------
                This file is named 'root_colors.inp' and is a text file that specifies all the input parameters to run color-color_plotter.py. The format of this file should be like root.inp as so:
                #all comments are preceded by a hashtag.
                #keywords will be followed by an '=' and then the value of the parameter. For example, for the directory name:
                dirname = dirname whatever it is
                #flag for whether or not to use special obj type
                spfg_[root] = 1 or 0
                coordinates = root_coords.csv
                #file formats to save plots as
                .pdf = y
                .jpg = n
                .png = n
                .eps = n
                #number of kmeans clusters
                nclusters = # of clusters
                #azim = ##, azimuthal angle to save 3d plot
                #elev = ##, elevation angle to save 3d plot
                #ccax = X-Y,Z-P, X-Y is X mag - Y mag on x axis, Z-P is Z mag - P mag on y axis
                #azim, elev, and ccax can be repeated with different values as many times as wanted
                
    ------- 2.7: Querying Images from IRSA -------
    * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/postage_stamp.py *
      root: root
      - Program Purpose: This program acquires images from IRSA. We are using GAIA EDR3 coordinates.
      - This program takes parameters specified in root-img.inp to get the images and coordinates from ohir_coords.csv.
      - Running Procedure: Create a file called root-img.inp in the same directory as the program. This .inp file is formatted as written below (2.7.1). Make sure you also have the root_coords.csv
      - Run this program from the same directory as get-dict-phot.py. This will be run from the terminal as: python postage_stamp.py root-img.inp

       
               ------- 2.7.1: Input Image Parameters File -------
                This file is named 'root-img.inp' and is a text file that specifies all the input parameters to run postage_stamp.py. The format of this file should be like root.inp as so:
                #all comments are preceded by a hashtag.
                #keywords will be followed by an '=' and then the value of the parameter. For example, for the directory name:
                dirname = dirname whatever it is
                
                #survey = Specifies the survey dataset(s) to be retrieved. ("IRIS" can be used instead of "IRAS" for the same image set.) It can be one dataset (e.g. "survey=SDSS") or multiple datasets (e.g. "survey=SDSS,DSS,SEIP"). It defaults to the five datasets (DSS, SDSS, 2MASS, WISE, IRAS). Here are all the catalogs you can use: DSS, SDSS, 2MASS, WISE, SEIP, AKARI, IRAS. IRAS is not very informative at this spatial scale and will be removed from the default in the future. VERY IMPORTANT: MAKE SURE THERE ARE NO SPACES IN BETWEEN THE CATALOGS, ONLY COMMAS!

                #subsetsize = Specifies the cutout size of the retrieved images in arcmin. It should be between 0.1 and 60.0 arcmin.
                
-- old version of root-img.inp (v0) --                
                #Here are the required parameters:
                #mission = mission name (see link above to get name of missions)
                #min_size = The minimum allowed cutout size for mission data in arcseconds.
                #max_size = The maximum allowed cutout size for mission data in arcseconds.
                #ntable_cutouts = The number of metadata tables to search, for cutouts of mission. Names of all N tables are listed below.
                #cutouttbln = Name of nth cutout table, dimensions of images (to get names, see the link above)
                #units = The units of the sizeX parameter
                #c_format = s (sexagesimal) or d (decimal)
                #sep = separator in your coords file
                #sizex = The image cutouts box size on the sky (units of this size parameter are specified by the next parameter, called "units", which can be deg, arcmin or arcsec). The size can be any number larger than zero (interpreted as a double) and smaller than the size specified by the "max_size" parameter (note units may be different). Note, in most cases, the maximum allowed sizeX value is 2.0 degrees, but it varies for the different data collections.
                
               ------- 2.7.2 - Parameter Records (from older version of root-img.inp (v0)) ----------------
               mission = GLIMPSE
               min_size = 1
               max_size = 599
               ntable_cutouts = 8
               cutouttbl1 = IRAC1_0.6ResMos
               cutouttbl2 = IRAC2_0.6ResMos
               cutouttbl3 = IRAC3_0.6ResMos
               cutouttbl4 = IRAC4_0.6ResMos
               cutouttbl5 = IRAC1_1.2ResMos
               cutouttbl6 = IRAC2_1.2ResMos
               cutouttbl7 = IRAC3_1.2ResMos
               cutouttbl8 = IRAC4_1.2ResMos
               units = arcsec
               c_format = s
               sep = ,
               sizex = 25
               
               
               mission = SAGE
               min_size = 756
               max_size = 7200
               ntable_cutouts = 3
               cutouttbl1 = mips24
               cutouttbl2 = mips70
               cutouttbl3 = mips160
               units = arcsec
               c_format = s
               sep = ,
               sizex = 900
               
               mission = LH
               min_size = 1
               max_size = 3600
               ntable_cutouts = 3
               cutouttbl1 = LH_J_band
               cutouttbl2 = LH_H_band
               cutouttbl3 = LH_K_band
               units = arcsec
               c_format = s
               sep = ,
               sizex = 240
               
               mission = SGOODS
               min_size = 1
               max_size = 600
               ntable_cutouts = 4
               cutouttbl1 = Spitzer_IRAC_DR3
               cutouttbl2 = Spitzer_MIPS_DR3
               cutouttbl3 = Spitzer_IRAC_DR2
               cutouttbl4 = Ancillary_Opical
               units = arcsec
               c_format = s
               sep = ,
               sizex = 15
               
               mission = GLIMPSE
               min_size = 1
               max_size = 599
               ntable_cutouts = 8
               cutouttbl1 = IRAC1_0.6ResMos
               cutouttbl2 = IRAC2_0.6ResMos
               cutouttbl3 = IRAC3_0.6ResMos
               cutouttbl4 = IRAC4_0.6ResMos
               cutouttbl5 = IRAC1_1.2ResMos
               cutouttbl6 = IRAC2_1.2ResMos
               cutouttbl7 = IRAC3_1.2ResMos
               cutouttbl8 = IRAC4_1.2ResMos
               units = arcsec
               c_format = s
               sep = ,
               sizex = 300

    ------- 2.8: Making Postage Stamp Mosaics -------
    * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/mosaic.py *
      root: root
      - Program Purpose: This program takes the outputs from the postage_stamp.py program and creates a mosaic of a user-specified size (specify in root-mosaic.inp) as well as a .csv for each star containing the default min and max intensities (if usr_flg = 0 (see .inp file)). We are using GAIA EDR3 coordinates.
      - This program takes parameters specified in root-mosaic.inp to get the images and coordinates from root_coords.csv. If usr_flg = 1, the program will also take the min, max intensities, and stretch specified in the .csv to edit the noise in the images.
      - Running Procedure: Create a file called root-mosaic.inp in the same directory as the program. This .inp file is formatted as written below (2.8.1). Make sure you also have the root_coords.csv
      - Run this program from the same directory as get-dict-phot.py. This will be run from the terminal as: python mosaic.py root-mosaic.inp
      
               ------- 2.8.1: Input Mosaic Parameters File -------
                This file is named 'root-mosaic.inp' and is a text file that specifies all the input parameters to run mosaic.py. The format of this file should be like root.inp as so:
                #all comments are preceded by a hashtag.
                #keywords will be followed by an '=' and then the value of the parameter. For example, for the directory name:
                dirname = dirname whatever it is

                #length = the length (x-axis) you want each of your images to be in arcseconds
                #width = the width (y-axis) you want each of your images to be in arcseconds
                #circle = the size of the circle that encircles the Gaia EDR3 coordinates in the image
                #c_format = s (sexagesimal) or d (decimal) of your coordinates in the coords file
                #usr_flg = when adjusting the min, max, and stretch of the images, do you want to use the default min, max, stretch (0) or your own (1)?
                #noplot = surveys that you DON'T want to plot (write each survey in capital letters and separate with a comma). If you would like to plot all of them, write "none"
                
    ------- 2.9: Post-RGB Star Candidate List -------
    * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/prgb_candidates.sh *
      root: root
      - Program Purpose: This program takes outputs from phot_plotter.py (specifically the dist_lum.csv) and get-dict-phot.py (specifically the phot.csv) and creates one giant list of post-RGB star candidates by seeing if their actual luminosities match the given criteria: >= 30 solar luminosities but <= 1500 solar luminosities. In addition, the OH maser line velocities are added to the list (which were found when get-dict-phot.py was ran).
      - Run this program in the same directory as your outputs from phot_plotter.py and get-dict-phot.py. If this is the first time you are running the program, type chmod +x prgb_candidates.sh, hit enter, then ./prgb_candidates.sh root (ex: if the root name of my files are abc, then I would type ./prgb_candidates.sh abc). Subsequent times, you can just type ./prgb_candidates.sh root
      
    ------- 2.10 - Creating LaTeX Tables Containing SED Input Data -------
    * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/pht_mod.py
    - Program Purpose: This program takes the sed2_name_raw5.csv that was created when get-dict-phot.py was ran and makes it into a LaTeX table (name_pht_MoD.tex). Specifically, the table outputs the following columns in the following order: filter name, wavelength, flux, flux error, offset, and catalog of name
    - Run this program in the same directory as your outputs from get-dict-phot.py. In the terminal, type python pht_mod.py name. What I mean by name is that if you look at the sed2 file names for your star, you will notice that there are underscores that replace the spaces. Make sure you write down the exact way the star is written in the file names. Ex: If I want to create a table for the star IRAS 23361+6437, I would type python pht_mod.py IRAS_23361+6437 because that is how the sed2 file name is written in the sed-for-MoD folder in ohir20-2.5asec.
    
    ------- 2.11 - Creating LaTeX Tables Containing Bailer-Jones/GAIA Data     -------
    * Location: /mchang-JPL2021/python_phot-lum-colors/photometry/bailjons_motion_data.py
    - Program Purpose: This program goes into the file ohir-all-final_Thres_eDR3-bailjons-cln.csv and extracts the name of source, input RA, input dec, output ra, output dec, offset, parallax, parallax error, proper motion, RA proper motion, RA proper motion error, dec proper motion, dec proper motion error, distance, lower limit distance, and upper limit distance of select stars (bailjons_motion_stars.inp) and outputs it into a .tex file formatted as a landscape table (bailjons_table.tex).
    - Run this program wherever you want but I did it in the photometry folder as that's where the other .py programs are. Make sure you have the correct file paths for your files. Have your list of stars ready (bailjons_motion_stars.inp). To run this program, type python bailjons_motion_data.py
    
               ------- 2.11.1: Input Stars Bailer-Jones/GAIA File -------
               This file is called bailjons_motion_stars.inp. It contains the list stars that you wish to extract info from the ohir-all-final_Thres_eDR3-bailjons-cln.csv file. The format of the file should be like this:
               
               #Comments should be preceded by a hashtage
               IRAS 23361+6437
               IRAS 18027-1209
               ...
               
               VERY IMPORTANT: Make sure the names of the stars are EXACTLY how you find it in ohir-all-final_Thres_eDR3-bailjons-cln.csv (using correct cases, no extra spaces, etc.), so the program can find the star.
    
======== 3 - Modeling SEDs of post-RGB stars ========
Make sure you have read the procedures on how to download MoD and copy directories and rename files for each of you stars.

From your outputs from get-dict-phot.py, go to sed-for-MoD and extract the files that have sed2_name.dat and place them in the appropriate directory in MoD (ex: if I have a .dat named sed2_IRAS_00007+5524.dat, I would place it in the folder named iras00007+5524 in my MoD directory). Rename the name part of the sed2_name.dat so it matches the root of the directory and all the other files. Then, go to the file root_in. This file contains the parameters for modeling the SED. This information can also be found in the MoD manual section 5: Input to MoD. Your root_in looks like this:

oh003.234
0.3
../ModelsDust/dhs0.7_a0.15amc0.95sic0.05mgs0.0.dat
../ModelsDust/s5500_g+0.5_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.dusty
5500.
0.335
1
1600.
-1.0 
1.0  
13. 1
50.   1
75.   1
2.0    0
!!!!!!!!!!!!!!!!!
V16a    
Lum tau Tc p

star=oh003.234
unique=`date +%Y%j%H%M`
PlotMyDusty.exe   ${star}_in
cp ${star}.inp    ${star}.inp.${unique}
cp ${star}_sed.ps ${star}_sed.ps.${unique}
more  ${star}.inp

The parameters are located above the exclamation points. Line by line, they are:

oh003.234 - An arbitrary string, but the root name of additional input file names
0.3 - Interstellar reddening, AV
../ModelsDust/dhs0.7_a0.15amc0.95sic0.05mgs0.0.dat - Name of the file containing the absorption and scattering coefficients (following DUSTY convention)
../ModelsDust/s5500_g+0.5_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.dusty - Name of the file containg the stellar input spectrum (following DUSTY convention)
5500. - Effective temperature (K)
0.335 - Distance in kpc
1 - Number of shells, N
1600. - Outer radius (radii), y(N) of these shells
-1.0 - Exponent of the density law, p(N), of these shells. Note that p(1) is one of the fit parameters in the current version of the code, so its value can be arbitrary (and put to -1 here)
1.0 - scaling factors, s(N)
13. 1 - free parameter 1, luminosity (solar units) and if it is fixed (0), or fitted (1).
50.   1 - free parameter 2, optical depth at 0.55 micron and if it is fixed (0), or fitted (1). As the log value is the parameter that is minimised, do not put 0 here, but a small number if you want to to start the minimisation from a very small optical depth value.
75.   1 - free parameter 3, temperature at inner dust radius (Kelvin) and if it is fixed (0), or fitted (1).
2.0    0 - free parameter 4, slope of the density law, p(1), and if it is fixed (0), or fitted (1).

In addition, make sure the line star=name has the root name of your star.

With each test run, record the results!

Running Procedure: First make sure .exe ready by running one of the following commands depending if you have pgplot or not: 
gfortran  -fbounds-check -O -o MyDusty.exe  MyDusty_v16a_noPGP.f -std=legacy (no pgplot)
gfortran  -fbounds-check -O -o MyDusty.exe  MyDusty_v16a.f -lpgplot -std=legacy (pgplot)

When done, then type: ./MyDusty.exe root_in

When finished, then type: gnuplot -p PlotMyDusty_oh003.234.gnu (if you need gnuplot, brew install gnuplot)
