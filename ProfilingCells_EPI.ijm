dir=getDirectory("image");
n = roiManager("count");
for (k=0; k<n; k++) {
	roiManager("select", k);
    roiManager("rename", k);
Stack.setPosition(1, 1, 1)
run("Clear Results");
run("Plot Profile");
Plot.getValues(x, y); 
  run("Clear Results"); 
  for (i=0; i<x.length; i++) { 
     setResult("x", i, x[i]); 
     setResult("y", i, y[i]); 
  } 
  setOption("ShowRowNumbers", false); 
  updateResults; 
saveAs("Measurements", dir+"Profiles//" + "Gray" + k + ".txt");
close();
Stack.setPosition(2, 1, 1)
run("Clear Results");
run("Plot Profile");
Plot.getValues(x, y); 
  run("Clear Results"); 
  for (i=0; i<x.length; i++) { 
     setResult("x", i, x[i]); 
     setResult("y", i, y[i]); 
  } 
  setOption("ShowRowNumbers", false); 
  updateResults; 
saveAs("Measurements", dir+"Profiles//" + "Red" + k + ".txt");
close();
Stack.setPosition(3, 1, 1)
run("Clear Results");
run("Plot Profile");
Plot.getValues(x, y); 
  run("Clear Results"); 
  for (i=0; i<x.length; i++) { 
     setResult("x", i, x[i]); 
     setResult("y", i, y[i]); 
  } 
  setOption("ShowRowNumbers", false); 
  updateResults; 
saveAs("Measurements", dir+"Profiles//" + "Blue" + k + ".txt");
close();
}
