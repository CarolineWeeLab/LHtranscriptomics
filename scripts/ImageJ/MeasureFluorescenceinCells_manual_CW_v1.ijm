// First, open an image and use multipoint TOOL to select all points. 
// Remember to click on CENTER of each cell at correct Z! 
// Remember to change oval size if necessary
// Need to download latest FIJI version for Z-information to be accurate
// Output will be in "Results", ordered by each cell (and all channels per cell)

run("Clear Results");
run("Set Measurements...", "mean center stack nan redirect=None decimal=3");
run("Measure");

ovalsize = 20; //sets number of pixels of circle diameter)
npoints = nResults;
nchannels = 4;

xvalue = newArray(npoints);
yvalue = newArray(npoints);
zvalue = newArray(npoints);

// to get the x y z coordinates for all points selected

for (i=0; i<npoints; i++) {
    x = getResult("X",i);
    y = getResult("Y", i);
    z = getResult("Slice", i);
    xvalue[i]= x;
    yvalue[i]= y;
    zvalue[i] = z;
}

run("Clear Results"); //now is when I make the actual measurements 

// to draw ovals over each point, and measure the intensity values

for (i=0; i<npoints; i++) {
    Stack.setChannel(1);	
    Stack.setSlice(zvalue[i]) 
    toUnscaled(xvalue[i], yvalue[i]);
    makeOval(xvalue[i]-ovalsize/2, yvalue[i]-ovalsize/2, ovalsize, ovalsize);

	run("Measure");
	setColor("gray");
	fill();

    //measure other channels
	for (j=2; j<nchannels+1; j++) {
	Stack.setChannel(j);
	run("Measure");
	}
	
	run("Select None");
	Stack.setChannel(1);
}
