//drag and drop CFW/Katushka single cell stacks obtained from macro 'manual selection of single cells and Katushka mask to CFW signal.ijm'
//then activate this macro to analyse CFW ROIs on CFW channels
//this macro analyses CFW level of all cells
number = nImages;
print(number);

//splits CFW and Katushka channels
for (i=1;i<number+1;i++) { 
		print(i);
        selectWindow(i+".tif"); 
        run("Split Channels");
        }

//duplicates CFW stacks only
number = nImages/2+1;
print(number);
for (i=1;i<number;i++) { 
		print(i);
		title1 = "C1-"+i+".tif";
		selectWindow(title1);
		run("Duplicate...", "duplicate");
		rename("a"+title1);
}

//makes ROI on CFW channels and measures ROI on duplicated CFW channels in time
for (i=1;i<number;i++) { 
		print(i);
		title1 = "C1-"+i+".tif";
		quotedTitle = "'" + title1 + "'";
		print(title1);
        title2 = "aC1-"+i+".tif"; 
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