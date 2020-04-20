//first open hyperstack with 3 channels (CFW, Katushka, BF) using virtual stack
//this macro analyses CFW level of living cells only

//selects hyperstack and sums slices
imageTitle=getTitle();
run("Z Project...", "projection=[Sum Slices] all");
//closes virtual stack window
selectWindow(imageTitle);
close();

//duplicates CFW (1) and Katushka (2) channel
selectWindow("SUM_" + imageTitle);
run("Duplicate...", "duplicate channels=1-2");//removes Bright Field channel and selects specific time frame optional
selectWindow("SUM_" + imageTitle);
close();

//duplication of Katushka channel and sums slices in time to generate 'foot print' of cell
selectWindow("SUM_" + imageTitle + "-1");
run("Duplicate...", "duplicate channels=2");
selectWindow("SUM_"+imageTitle+"-2");//selects Katushka channel
run("Z Project...", "projection=[Sum Slices]");//sum slices in time of the Katushka channel
selectWindow("SUM_"+imageTitle+"-2");
close();

waitForUser("Mannually select single cells and add them to ROI manager by selection and pressing 't'");

selectWindow("SUM_SUM_"+imageTitle+"-2");
close();

//counts number of ROIs in ROI manager
nROIs = roiManager("count");
print("nROIs = ");
print(nROIs);

//duplicates time series of each individual cell
for (i = 0; i < nROIs; i++)
{
selectWindow("SUM_" + imageTitle + "-1");
roiManager("Select", i);
run("Duplicate...", "duplicate");
}

selectWindow("SUM_" + imageTitle + "-1");
close();

waitForUser("Mannually go through single cells and remove any double cells/debris from images");

// get image IDs of all open images and saves them to user specified folder
dir = getDirectory("Choose a Directory"); 
ids=newArray(nImages); 
for (i=0;i<nImages;i++) { 
        selectImage(i+1); 
        title = getTitle; 
        print(title); 
        ids[i]=getImageID; 
	
        saveAs("tiff", dir+i+1); // this saves 1 to infinity instead of 'getTitle'. Since macro requires short file names??
}

selectWindow("ROI Manager");
      run("Close");

number = nImages;
print(number);

//splits CFW and Katushka channels
for (i=1;i<number+1;i++) { 
		print(i);
        selectWindow(i+".tif"); 
        run("Split Channels");
        }

//makes ROI on Katushka channels and measures ROI on CFW channels in time
for (i=1;i<number+1;i++) { 
		print(i);
		title1 = "C1-"+i+".tif";
		quotedTitle = "'" + title1 + "'";
		print(title1);
        title2 = "C2-"+i+".tif";
        selectWindow(title2); 
		run("8-bit");
		setAutoThreshold("Default dark");
		//run("Threshold...");
		run("Convert to Mask", "method=Default background=Dark calculate");
		//run("Close");
		run("Set Measurements...", "area mean display redirect="+quotedTitle);
		run("Analyze Particles...", "size=2-Infinity display add stack");
		//updateResults();
}



//closes all open images
while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 

 	selectWindow("Log");
      run("Close");
      selectWindow("ROI Manager");
      run("Close");

waitForUser("Mannually copy paste results to excel");
selectWindow("Results");
      run("Close");
