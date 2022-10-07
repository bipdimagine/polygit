
var dataStorePlotValue;
var dataStoreGenesValue;

var gridPlot;
var tabPos ="";
var tabRatio = "";


//----------------------------------
// plot simple
//----------------------------------
function Plotsimple()
{
	var url_plot = url_path + "/manta/Plot.pl";
	
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	
	var filename =par0[1];
	dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plot+ "?ficin=" + filename });
	
	dataStorePlotValue.fetch({
			onComplete: function(items,request){
                	var item = items[0];
  
   					var listPos = dataStorePlotValue.getValue(item, "POS"); 
   					var listValues = dataStorePlotValue.getValue(item, "VALUES");
 					
 					var tabPos = listPos.split(" ");
   					var tabValues = listValues.split(" ");
 					
 		
  					var trace1 = {
  							x:tabPos,
  							y:tabValues,						
  							mode: 'markers',
  							showlegend:false,
  							name:' titi',
  							marker: {
      								color: 'grey',
      								size:3
    						}
					};
		
					data=[trace1];
					layout = {
  							xaxis : {title : 'position'},
  							yaxis:   {title : 'value'},
						};
						
				Plotly.newPlot('plot1', data, layout);
		}
  	});		
  	
}


/////////////////////////////////////////////////////////////////////////////////////
//      pour les plots des frequences alleliques seules
/////////////////////////////////////////////////////////////////////////////////////

function launch_plotFreqAllelique(chr)
{
	var out;
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
	
	var project =par0[1];
   	var patient = par1[1];
   	
   	if(chr == 23) {chr="X";}
   	if(chr == 24) {chr="Y";}
   	
   	var url_plot;
	
	url_plot = "PolyCytoPlotFreqAllelique.html?project=" + project + "&patient=" + patient + "&chr=" + chr ;
	var myWindow = window.open(url_plot,"_blank",""); 
}
 
 function PlotFreqAllelique()
{
	var parameters = location.search.split("&");
	
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");	 
	var par2 = parameters[2].split("=");
	
	var project =par0[1];
	var patient =par1[1];
	var chr=par2[1];
	
	var start = 0;
	var end = 0;
	
	if (parameters.length == 4)
	{
			var par3 = parameters[3].split("=");
			var pos =par3[1];
			var start = parseInt(pos) - 1000000;
			var end =  parseInt(pos) + 1000000;
	}
	
	// plot = la frequence allelique 
	PlotValuesFreqAllelique(project,patient,chr,start,end);
	setTimeout(function(){document.getElementById('labw1').style.visibility = 'hidden';},5000);
 }
 
 function PlotValuesFreqAllelique(project,patient,chr,start,end)
{
	var url_plot = url_path + "/manta/PolyCytoPlotFreqAll.pl";
	var titre;
	
	if (start == 0)
	{
		dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plot + "?project=" + project + "&patient=" + patient + "&chr=" + chr + "&type=" + 0});
		titre = '<br><label style="font-size:16px; color:orange">' +patient + ' / chr' + chr + '</label><br><br>';
	}
	else
	{
		dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plot + "?project=" + project + "&patient=" + patient + "&chr=" + chr + "&start=" + start + "&end=" + end + "&type=" + 0});
		titre = '<br><label style="font-size:16px; color:orange">' +patient + ' / chr' + chr + '_' + start + '_' + end + '</label><br><br>';
	}
	
	
	document.getElementById('idtitre').innerHTML = titre;
	
	var trio;
	
	var tabPosBA;
   	var tabMother;
   	var tabFather;
   	var tabBAX=[];
   	var tabBAYM=[];
   	var tabBAYF=[];
	var tabplotBAX=[start,end];
   	var tabplotBAYM=[0,100];
   	var tabplotBAYF=[0,100];
   	
   	var tabPosBA2;
   	var tabMother2;
   	var tabFather2;
   	var tabBAX2=[];
   	var tabBAYM2=[];
   	var tabBAYF2=[];
	var tabplotBAX2=[start,end];
   	var tabplotBAYM2=[0,100];
   	var tabplotBAYF2=[0,100];
   	
   	
	dataStorePlotValue.fetch({
			onComplete: function(items,request){
                	var item = items[0];
  
   					tabPosBA = dataStorePlotValue.getValue(item, "POSball"); 
   					tabMother = dataStorePlotValue.getValue(item, "MOTHER");
 					tabFather = dataStorePlotValue.getValue(item, "FATHER");
 					
 					tabPosBA2 = dataStorePlotValue.getValue(item, "POSball_2"); 
   					tabMother2 = dataStorePlotValue.getValue(item, "MOTHER_2");
 					tabFather2 = dataStorePlotValue.getValue(item, "FATHER_2");

   					var tabBAX = tabPosBA.split(" ");
   					var tabBAYM = tabMother.split(" ");
   					var tabBAYF = tabFather.split(" ");
   					
   					var tabBAX2 = tabPosBA2.split(" ");
   					var tabBAYM2 = tabMother2.split(" ");
   					var tabBAYF2 = tabFather2.split(" ");
   					
   						for ( var i = 0 ; i < tabBAX.length ; i++ ) {
  									tabplotBAX.push(tabBAX[i]);
  									tabplotBAYM.push(tabBAYM[i]);			
  									tabplotBAYF.push(tabBAYF[i]);
  						}
	
  						var trace1 = {
  							x:tabplotBAX,
  							y:tabplotBAYM,						
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'pink',
      								size:4
    						}
						};
				
						var trace2 = {
  							x:tabplotBAX,
  							y:tabplotBAYF,							
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'lightblue',
      								size:4
    						}
						};
						
						for ( var i = 0 ; i < tabBAX2.length ; i++ ) {
  									tabplotBAX2.push(tabBAX2[i]);
  									tabplotBAYM2.push(tabBAYM2[i]);			
  									tabplotBAYF2.push(tabBAYF2[i]);
  						}
	
  						var trace3 = {
  							x:tabplotBAX2,
  							y:tabplotBAYM2,						
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'deeppink',
      								size:4
    						}
						};
				
						var trace4 = {
  							x:tabplotBAX2,
  							y:tabplotBAYF2,							
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
						};
						
						
						data=[trace3,trace4,trace1,trace2];
						layout = {
  								xaxis : {title : 'genomic position'},
  								yaxis:   {domain: [0.4, 1], title : 'Allelic balance', tickmode: "array", tickvals: [0, 25, 50, 75, 100], ticktext: ['0', '25', '50', '75', '100']},
						};
						
				var config = { modeBarButtonsToRemove: ['pan2d', 'select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d','toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'], displayModeBar: true};
				Plotly.newPlot('plot1', data, layout,config);
		}
  	});		
  	
}


////////////////////////////////////////////////////////
// pour les plots des chromosomes
/////////////////////////////////////////////////////////

function launch_plotChr(chr)
{
	var out;
	
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
	
	var project =par0[1];
   	var patient = par1[1];
   	
   	if(chr == 23) {chr="X";}
   	if(chr == 24) {chr="Y";}
   	
   	var url_plot;
	
	url_plot = "PolyCytoPlotChromosome.html?project=" + project + "&patient=" + patient + "&sens=H" + "&chr=" + chr ;
	var myWindow = window.open(url_plot,"_blank",""); 
}

function launch_plotAllChr()
{
	var out;
	
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
	
	var project =par0[1];
   	var patient = par1[1];
   	
   	var url_plot;
   
	url_plot = "PolyCytoPlotAllChromosome.html?project=" + project + "&patient=" + patient;
	var myWindow = window.open(url_plot,"_blank",""); 
}


 function PlotChr()
{
	var parameters = location.search.split("&");
	
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");	 
	var par2 = parameters[2].split("=");
	var par3 = parameters[3].split("=");
	
	var project =par0[1];
	var patient =par1[1];
	var sens = par2[1];
	var chr=par3[1];
	
	// plot = les ratios wisecondor 
	if (sens == "V")
	{
		PlotValuesChrVert(project,patient,chr);
	}
	else
	{
		document.getElementById('divbody').innerHTML="<br>";
		PlotValuesChr(project,patient,chr);
		document.getElementById('divbody').innerHTML="<br>";
	}
	setTimeout(function(){document.getElementById('labw1').style.visibility = 'hidden';},5000);
 }
 
 function PlotValuesChr(project,patient,chr)
{
	var url_plotChr = url_path + "/manta/PolyCytoPlotChromosome.pl";
	dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plotChr + "?project=" + project + "&patient=" + patient + "&chr=" + chr });
	
	var trio;
	var type=0;
	
	var tabPosR;
	var tabRatio;
	var tabZScore;
	var tabRX;
	var tabRY;
	var tabRZ;
	
	var tabPosR_dad;
	var tabRatio_dad;
	var tabZScore_dad;
	var tabRX_dad;
	var tabRY_dad;
	var tabRZ_dad;
	
	var tabPosR_mom;
	var tabRatio_mom;
	var tabZScore_mom;
	var tabRX_mom;
	var tabRY_mom;
	var tabRZ_mom;
	
	var tabplotRX=[0,0];
   	var tabdelRX=[];
   	var tabdupRX=[];
   	var tabplotRY=[-1,1];
   	var tabdelRY=[];
   	var tabdupRY=[];
   	
   	
   	var tabplotRX_dad=[0,0];
   	var tabdelRX_dad=[];
   	var tabdupRX_dad=[];
   	var tabplotRY_dad=[-1,1];
   	var tabdelRY_dad=[];
   	var tabdupRY_dad=[];
   	
   	var tabplotRX_mom=[0,0];
   	var tabdelRX_mom=[];
   	var tabdupRX_mom=[];
   	var tabplotRY_mom=[-1,1];
   	var tabdelRY_mom=[];
   	var tabdupRY_mom=[];
   	
   	
   var graphDiv = document.getElementById('plotdiv');	
   	
	dataStorePlotValue.fetch({
			onComplete: function(items,request){
                	var item = items[0];
  					
  					type = parseInt(dataStorePlotValue.getValue(item, "TYPE"));		
   					trio = parseInt(dataStorePlotValue.getValue(item, "PLOT"));			
   					tabPosR = dataStorePlotValue.getValue(item, "POSratio"); 
   					tabRatio = dataStorePlotValue.getValue(item, "RATIO");
  					tabZScore = dataStorePlotValue.getValue(item, "Zscore");
  					
   					tabRX = tabPosR.split(" ");
   					tabRY = tabRatio.split(" ");
   					tabRZ = tabZScore.split(" ");
   					
   					for ( var i = 0 ; i < tabRX.length ; i++ ) 
   					{
   						// points aberrants
   						if ( tabRY[i]  < -2 ) 
   						{
 							tabRY[i] = -2;
 						}
 						
 						if ( tabRY[i]  > 2) 
   						{
 							tabRY[i] = 2;
 						}
 					
 						if ( ( tabRZ[i] < -5 ) || ( tabRZ[i] > 5) )
 						{
  							if  ( tabRY[i] >0.22 ) 
  							{
  								tabdupRX.push(tabRX[i]);
  								tabdupRY.push(tabRY[i]);
  							}
  							else
  							{
  								if ( tabRY[i]  < -0.26 ) 
  								{
  									tabdelRX.push(tabRX[i]);
  									tabdelRY.push(tabRY[i]);
  								}
  							}
  						}
  						else
  						{
  							tabplotRX.push(tabRX[i]);
  							tabplotRY.push(tabRY[i]);
  						}
  					}
  					
  					var trace1 = {
  							x:tabplotRX,
  							y:tabplotRY,
  							yaxis: 'y1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
					};
					
						
					var trace2 = {
  							x:tabdupRX,
  							y:tabdupRY,
  							yaxis: 'y1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
					};
					
					
					var trace3 = {
  							x:tabdelRX,
  							y:tabdelRY,
  							yaxis: 'y1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'red',
      								size:4
    						}
					};
	
	
					data=[trace1,trace2,trace3];
					layout = {
  							yaxis1:   {title : 'Wc ratio : Child'},
  							xaxis : 	  {title : 'Chr' + chr},
					};
					
					if (trio)
					{
						// ratio des parents
						
						// 1 pour le pere
						tabPosR_dad = dataStorePlotValue.getValue(item, "POSratio_dad"); 
   					 	tabRatio_dad = dataStorePlotValue.getValue(item, "RATIO_dad");
   					 	tabZScore_dad = dataStorePlotValue.getValue(item, "Zscore_dad");
  	
   						tabRX_dad = tabPosR_dad.split(" ");
   						tabRY_dad= tabRatio_dad.split(" ");
 						tabRZ_dad = tabZScore_dad.split(" ");
 						
 						for ( var i = 0 ; i < tabRX_dad.length ; i++ ) 
   						{
   							//points aberrants
   							if ( tabRY_dad[i]  < -2 ) 
   							{
 								tabRY_dad[i] = -2;
 							}
 						
 							if ( tabRY_dad[i]  > 2) 
   							{
 								tabRY_dad[i] = 2;
 							}
   						
   							if ( ( tabRZ_dad[i] < -5 ) || ( tabRZ_dad[i] > 5) )
 							{
   								if  ( tabRY_dad[i] >0.22 ) 
  								{
  									tabdupRX_dad.push(tabRX_dad[i]);
  									tabdupRY_dad.push(tabRY_dad[i]);
  								}
  								else
  								{
  									if ( tabRY_dad[i]  < -0.26 ) 
  									{
  										tabdelRX_dad.push(tabRX_dad[i]);
  										tabdelRY_dad.push(tabRY_dad[i]);
  									}
  								}
  							}
  							else
  							{
  								tabplotRX_dad.push(tabRX_dad[i]);
  								tabplotRY_dad.push(tabRY_dad[i]);
  							}
  						}

  						var trace4 = {
  							x:tabplotRX_dad,
  							y:tabplotRY_dad,
  							yaxis: 'y2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
						};
	
						
						var trace5 = {
  							x:tabdupRX_dad,
  							y:tabdupRY_dad,
  							yaxis: 'y2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
						};
						
						var trace6 = {
  							x:tabdelRX_dad,
  							y:tabdelRY_dad,
  							yaxis: 'y2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'red',
      								size:4
    						}
						};
						
						// pour la mere
						tabPosR_mom = dataStorePlotValue.getValue(item, "POSratio_mom"); 
   					 	tabRatio_mom = dataStorePlotValue.getValue(item, "RATIO_mom");
   					 	tabZScore_mom = dataStorePlotValue.getValue(item, "Zscore_mom");
  
   						tabRX_mom = tabPosR_mom.split(" ");
   						tabRY_mom= tabRatio_mom.split(" ");
   						tabRZ_mom = tabZScore_mom.split(" ");
   						
 						for ( var i = 0 ; i < tabRX_mom.length ; i++ ) 
   						{
   							
   							//points aberrants
   							if ( tabRY_mom[i]  < -2 ) 
   							{
 								tabRY_mom[i] = -2;
 							}
 						
 							if ( tabRY_mom[i]  > 2) 
   							{
 								tabRY_mom[i] = 2;
 							}
   						
   							if ( ( tabRZ_mom[i] < -5 ) || ( tabRZ_mom[i] > 5) )
 							{
   								if  ( tabRY_mom[i] >0.22 ) 
  								{
  									tabdupRX_mom.push(tabRX_mom[i]);
  									tabdupRY_mom.push(tabRY_mom[i]);
  								}
  								else
  								{
  									if ( tabRY_mom[i]  < -0.26 ) 
  									{
  										tabdelRX_mom.push(tabRX_mom[i]);
  										tabdelRY_mom.push(tabRY_mom[i]);
  									}
  								}
  							}
  							else
  							{
  								tabplotRX_mom.push(tabRX_mom[i]);
  								tabplotRY_mom.push(tabRY_mom[i]);
  							}
  						}
 
  						var trace7 = {
  							x:tabplotRX_mom,
  							y:tabplotRY_mom,
  							yaxis: 'y3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
						};
						
						var trace8 = {
  							x:tabdupRX_mom,
  							y:tabdupRY_mom,
  							yaxis: 'y3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
						};
						
						var trace9 = {
  							x:tabdelRX_mom,
  							y:tabdelRY_mom,
  							yaxis: 'y3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'red',
      								size:4
    						}
						};
						
						data=[trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9];
						
						layout = {
  								yaxis1:   {domain: [0.0, 0.3], title : 'EA'},
  								yaxis2:   {domain: [0.35, 0.65], title : ' F', color:'blue'},
  								yaxis3:   {domain: [0.7, 1], title : 'M', color:'deeppink'},
  								xaxis : 	  {title : 'Chr' + chr},
						};
				} // fin if trio
	
				var config  = { modeBarButtonsToRemove: ['pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d','toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'], displayModeBar: true};
				
				if (type == 0) {
  						document.getElementById("plotdiv").style="height:92%;width:98%; border: 4px solid lightgrey;";
  				}
  				else{
  						if (type == 1) {
  							document.getElementById("plotdiv").style="height:92%;width:98%; border: 4px solid blue;";
  						}
  						else{
  									if (type == 2){
  										document.getElementById("plotdiv").style="height:92%;width:98%; border: 4px solid red;";
  									}
  									else{
  										document.getElementById("plotdiv").style="height:92%;width:98%; border: 4px solid purple;"
  									}
  						}
  				}
  				
				Plotly.newPlot(graphDiv, data, layout,config);
				
				graphDiv.on('plotly_selected', function(eventData)
 				{
 
					var max = 0;
					var min = 250000000;
					
 					eventData.points.forEach(function(pt) {
  						min = Math.min(min,pt.x);
  						max = Math.max(max,pt.x);
 					});
 
 					var locus = chr+":"+min+"-"+max;
 					viewGenes(locus,project,patient);
 					
 					Plotly.restyle(graphDiv, data, layout,config);
				});
		}
  	});		
}
 
 
 
 function viewGenes(locus,project,patient)
 {
 	var username = dojo.cookie("username");
 	var url_view_gene = "/polyweb/vector/get_variants_from_all_my_projects.html?input=" + locus + "&project="+ project + "&patient="+patient + "&user="+username;
	var myWindow = window.open(url_view_gene,"_blank",""); 
 }
 
 function viewGeneByName(name,project,patient)
 {
 	//var url_view_gene = "PolyCytoTest.html?genes=" + name;
 	
 	var username = dojo.cookie("username");
 	var url_view_gene = "/polyweb/vector/get_variants_from_all_my_projects.html?input=" + name + "&project="+ project + "&patient="+patient + "&user="+username;
 	
	var myWindow = window.open(url_view_gene,"_blank",""); 
 }
 
 
 function getGeneInfos() 
{
	var username = dojo.cookie("username");
   	
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
   	
   	var what = par0[0];
 	var url_genes;
 	

	if (what == "?locus")
	{ 
   		var locus =par0[1];
   		url_genes = "/cgi-bin/polymorphism-cgi/json_output_nodb/check_genes_transcripts_locus.pl?user=username&locus="+locus;
   		
   		// recuperation du json
   		$.getJSON( url_genes, function( data ) {
  				
  				var myPanel = $(data.html);
				var myNode = document.getElementById("divgenes");
				
				myPanel.appendTo(myNode);
		});
   	}
   	
   	if (what == "?genes")
	{ 
   		var gene =par0[1];
   		url_genes = "/cgi-bin/polymorphism-cgi/json_output_nodb/check_genes_transcripts_locus.pl?user=username&genes="+gene;
   		$.getJSON( url_genes, function( data ) {
  				
  				var myPanel = $(data.html);
				var myNode = document.getElementById("divgenes");
				
				myPanel.appendTo(myNode);
		});
 	}
}


///////////////////////////////////////////////////////////////////
//   pour les plots chromosome a la vertical
///////////////////////////////////////////////////////////////////
 
 
 function PlotValuesChrVert(project,patient,chr)
{
	var url_plotChr = url_path + "/manta/PolyCytoPlotChromosome.pl";
	dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plotChr + "?project=" + project + "&patient=" + patient + "&chr=" + chr });
	
	var trio;
	var type=0;
	
	var tabPosR;
	var tabRatio;
	var tabZScore;
	var tabRX;
	var tabRY;
	var tabRZ;
	
	var tabPosR_dad;
	var tabRatio_dad;
	var tabZScore_dad;
	var tabRX_dad;
	var tabRY_dad;
	var tabRZ_dad;
	
	var tabPosR_mom;
	var tabRatio_mom;
	var tabZScore_mom;
	var tabRX_mom;
	var tabRY_mom;
	var tabRZ_mom;
	
	var tabplotRX=[0,0];
   	var tabdelRX=[];
   	var tabdupRX=[];
   	var tabplotRY=[-1,1];
   	var tabdelRY=[];
   	var tabdupRY=[];
   	
   	
   	var tabplotRX_dad=[0,0];
   	var tabdelRX_dad=[];
   	var tabdupRX_dad=[];
   	var tabplotRY_dad=[-1,1];
   	var tabdelRY_dad=[];
   	var tabdupRY_dad=[];
   	
   	var tabplotRX_mom=[0,0];
   	var tabdelRX_mom=[];
   	var tabdupRX_mom=[];
   	var tabplotRY_mom=[-1,1];
   	var tabdelRY_mom=[];
   	var tabdupRY_mom=[];
   	
   var Xmax=0;
   
   var graphDiv = document.getElementById('plotdiv');	
   	
	dataStorePlotValue.fetch({
			onComplete: function(items,request){
                	var item = items[0];
  					
  					type = parseInt(dataStorePlotValue.getValue(item, "TYPE"));		
   					trio = parseInt(dataStorePlotValue.getValue(item, "PLOT"));			
   					tabPosR = dataStorePlotValue.getValue(item, "POSratio"); 
   					tabRatio = dataStorePlotValue.getValue(item, "RATIO");
  					tabZScore = dataStorePlotValue.getValue(item, "Zscore");
  					
   					tabRX = tabPosR.split(" ");
   					tabRY = tabRatio.split(" ");
   					tabRZ = tabZScore.split(" ");
   					
   					
   					for ( var i = 0 ; i < tabRX.length ; i++ ) 
   					{
   						if (tabRX[i] > Xmax) 
   						{
   							Xmax = parseInt(tabRX[i]);
   						}
   						
   						// points aberrants
   						if ( tabRY[i]  < -2 ) 
   						{
 							tabRY[i] = -2;
 						}
 						
 						if ( tabRY[i]  > 2) 
   						{
 							tabRY[i] = 2;
 						}
 					
 						if ( ( tabRZ[i] < -5 ) || ( tabRZ[i] > 5) )
 						{
  							if  ( tabRY[i] >0.22 ) 
  							{
  								tabdupRX.push(tabRX[i]);
  								tabdupRY.push(tabRY[i]);
  							}
  							else
  							{
  								if ( tabRY[i]  < -0.26 ) 
  								{
  									tabdelRX.push(tabRX[i]);
  									tabdelRY.push(tabRY[i]);
  								}
  							}
  						}
  						else
  						{
  							tabplotRX.push(tabRX[i]);
  							tabplotRY.push(tabRY[i]);
  						}
  					}
  					
  					var trace1 = {
  							y:tabplotRX,
  							x:tabplotRY,
  							xaxis1: 'x1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
					};
					
						
					var trace2 = {
  							y:tabdupRX,
  							x:tabdupRY,
							xaxis1: 'x1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
					};
					
					
					var trace3 = {
  							y:tabdelRX,
  							x:tabdelRY,
  							xaxis1: 'x1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'red',
      								size:4
    						}
					};
	
	
					data=[trace1,trace2,trace3];
					layout = {
  							xaxis1:   {title : 'Wc ratio : Child'},
  							yaxis : 	  {title : 'Chr' + chr, range: [Xmax,0]},
					};
					
					if (trio)
					{
						// ratio des parents
						
						// 1 pour le pere
						tabPosR_dad = dataStorePlotValue.getValue(item, "POSratio_dad"); 
   					 	tabRatio_dad = dataStorePlotValue.getValue(item, "RATIO_dad");
   					 	tabZScore_dad = dataStorePlotValue.getValue(item, "Zscore_dad");
  	
   						tabRX_dad = tabPosR_dad.split(" ");
   						tabRY_dad= tabRatio_dad.split(" ");
 						tabRZ_dad = tabZScore_dad.split(" ");
 						
 						for ( var i = 0 ; i < tabRX_dad.length ; i++ ) 
   						{
   							//points aberrants
   							if ( tabRY_dad[i]  < -2 ) 
   							{
 								tabRY_dad[i] = -2;
 							}
 						
 							if ( tabRY_dad[i]  > 2) 
   							{
 								tabRY_dad[i] = 2;
 							}
   						
   							if ( ( tabRZ_dad[i] < -5 ) || ( tabRZ_dad[i] > 5) )
 							{
   								if  ( tabRY_dad[i] >0.22 ) 
  								{
  									tabdupRX_dad.push(tabRX_dad[i]);
  									tabdupRY_dad.push(tabRY_dad[i]);
  								}
  								else
  								{
  									if ( tabRY_dad[i]  < -0.26 ) 
  									{
  										tabdelRX_dad.push(tabRX_dad[i]);
  										tabdelRY_dad.push(tabRY_dad[i]);
  									}
  								}
  							}
  							else
  							{
  								tabplotRX_dad.push(tabRX_dad[i]);
  								tabplotRY_dad.push(tabRY_dad[i]);
  							}
  						}

  						var trace4 = {
  							y:tabplotRX_dad,
  							x:tabplotRY_dad,
  							xaxis: 'x2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
						};
	
						
						var trace5 = {
  							y:tabdupRX_dad,
  							x:tabdupRY_dad,
  							xaxis: 'x2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
						};
						
						var trace6 = {
  							y:tabdelRX_dad,
  							x:tabdelRY_dad,
  							xaxis: 'x2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'red',
      								size:4
    						}
						};
						
						// pour la mere
						tabPosR_mom = dataStorePlotValue.getValue(item, "POSratio_mom"); 
   					 	tabRatio_mom = dataStorePlotValue.getValue(item, "RATIO_mom");
   					 	tabZScore_mom = dataStorePlotValue.getValue(item, "Zscore_mom");
  
   						tabRX_mom = tabPosR_mom.split(" ");
   						tabRY_mom= tabRatio_mom.split(" ");
   						tabRZ_mom = tabZScore_mom.split(" ");
   						
 						for ( var i = 0 ; i < tabRX_mom.length ; i++ ) 
   						{
   							
   							//points aberrants
   							if ( tabRY_mom[i]  < -2 ) 
   							{
 								tabRY_mom[i] = -2;
 							}
 						
 							if ( tabRY_mom[i]  > 2) 
   							{
 								tabRY_mom[i] = 2;
 							}
   						
   							if ( ( tabRZ_mom[i] < -5 ) || ( tabRZ_mom[i] > 5) )
 							{
   								if  ( tabRY_mom[i] >0.22 ) 
  								{
  									tabdupRX_mom.push(tabRX_mom[i]);
  									tabdupRY_mom.push(tabRY_mom[i]);
  								}
  								else
  								{
  									if ( tabRY_mom[i]  < -0.26 ) 
  									{
  										tabdelRX_mom.push(tabRX_mom[i]);
  										tabdelRY_mom.push(tabRY_mom[i]);
  									}
  								}
  							}
  							else
  							{
  								tabplotRX_mom.push(tabRX_mom[i]);
  								tabplotRY_mom.push(tabRY_mom[i]);
  							}
  						}
 
  						var trace7 = {
  							y:tabplotRX_mom,
  							x:tabplotRY_mom,
  							xaxis: 'x3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
						};
						
						var trace8 = {
  							y:tabdupRX_mom,
  							x:tabdupRY_mom,
  							xaxis: 'x3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
						};
						
						var trace9 = {
  							y:tabdelRX_mom,
  							x:tabdelRY_mom,
  							xaxis: 'x3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'red',
      								size:4
    						}
						};
						
						data=[trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9];
						
						layout = {
  								xaxis1:   {domain: [0.0, 0.3], title : 'EA'},
  								xaxis2:   {domain: [0.35, 0.65], title : ' F', color:'blue'},
  								xaxis3:   {domain: [0.7, 1], title : 'M', color:'deeppink'},
  								yaxis : 	  {title : 'Chr' + chr, range: [Xmax,0]},
						};
				} // fin if trio
	
				var config = { modeBarButtonsToRemove: ['pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d','toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'], displayModeBar: true};
				
				if (type == 0) {
  						document.getElementById("plotdiv").style="height:98%;width:98%; border: 4px solid grey;";
  				}
  				else{
  						if (type == 1) {
  							document.getElementById("plotdiv").style="height:98%;width:98%; border: 4px solid blue;";
  						}
  						else{
  									if (type == 2){
  										document.getElementById("plotdiv").style="height:98%;width:98%; border: 4px solid red;";
  									}
  									else{
  										document.getElementById("plotdiv").style="height:98%;width:98%; border: 4px solid purple;"
  									}
  						}
  				}
  				
				Plotly.newPlot(graphDiv, data, layout,config);
				
				graphDiv.on('plotly_selected', function(eventData)
 				{
 
					var max = 0;
					var min = 250000000;
					
 					eventData.points.forEach(function(pt) {
  						min = Math.min(min,pt.y);
  						max = Math.max(max,pt.y);
 					});
 
 					var locus = chr+":"+min+"-"+max;
 					viewGenes(locus,project,patient);
 					
 					Plotly.restyle(graphDiv, data, layout,config);
				});
		}
  	});		
}

////////////////////////////////////// 
// pour les plots des CNV 
////////////////////////////////////////
function launch_plot(transmission,type,chr,debcnv,fincnv)
{
	var out;
	var parameters = location.search.split("&");
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");
	
	var project =par0[1];
   	var patient = par1[1];
   	
   	if(chr == 23) {chr="X";}
   	if(chr == 24) {chr="Y";}
   
	var url_plot;
	
	url_plot = "PolyCytoCNVPlot.html?project=" + project + "&patient=" + patient + "&chr=" + chr + "&debcnv=" + debcnv + "&fincnv=" + fincnv + "&type=" + type + "&transmission=" + transmission;
	var myWindow = window.open(url_plot,"_blank",""); 
}

 
 function Plots()
{
	var parameters = location.search.split("&");
	
	var par0 = parameters[0].split("=");
	var par1 = parameters[1].split("=");	 
	var par2 = parameters[2].split("=");
	var par3 = parameters[3].split("=");	 
	var par4 = parameters[4].split("=");
	var par5 = parameters[5].split("=");	
	var par6 = parameters[6].split("=");	 
	
	var project =par0[1];
	var patient =par1[1];
	var chr=par2[1];
	var debcnv = par3[1];
	var fincnv = par4[1];
	var type = par5[1];
	var transmission = par6[1];
	
	// plot = les ratios wisecondor 
	PlotValues(transmission,type,project,patient,chr,debcnv,fincnv);
	setTimeout(function(){document.getElementById('labw1').style.visibility = 'hidden';},5000);
 }
 	
function PlotValues(transmission,type,project,patient,chr,debcnv,fincnv)
{
	var url_plotCNV = url_path + "/manta/PolyCytoPlotCNV.pl";
	dataStorePlotValue = new dojo.data.ItemFileReadStore({ url: url_plotCNV + "?project=" + project + "&patient=" + patient + "&chr=" + chr + "&type=" + type + "&debcnv=" + debcnv  + "&fincnv=" + fincnv});

	var ttype; // type version text
	
	if (type == 1) { ttype = "Deletion";}
	if (type == 2) { ttype = "Duplication";}
	if (type == 3) { ttype = "Triplication";}
	if (type == 4) { ttype = "Quadruplication";}
	
	var colortype = "blue";
	if (ttype == "Deletion") {colortype = "red";}
	
	var colortrans = "red";
	if (transmission == "mother") {colortrans = "pink";}
	if (transmission == "father") {colortrans = "lightblue";}
	if (transmission == "both") {colortrans = "grey";}
	
	if (transmission =="-")
	{
		var titre = '<br><label style="font-size:16px; color:orange">' +patient + ' / chr' + chr + ':' + debcnv + '_' + fincnv +'</label><br><b><label style="font-size:16px; color: ' + colortype +'">' + ttype  + '</label></b><br><br>';
	}
	else
	{
		if(transmission == "denovo")
		{
			var titre = '<br><label style="font-size:16px; color:orange">' +patient + ' / chr' + chr + ':' + debcnv + '_' + fincnv +'</label><br><br><b><label style="font-size:16px; color: ' + colortype +'">' + ttype  + '</label>' + '<label style="font-size:16px; color:' + colortrans+'">&nbsp&nbsp' + transmission + '</label></b><br><br>';
		}
		else
		{
			var titre = '<br><label style="font-size:16px; color:orange">' +patient + ' / chr' + chr + ':' + debcnv + '_' + fincnv +'</label><br><br><b><label style="font-size:16px; color: ' + colortype +'">' + ttype  + '</label>' + '<label style="font-size:16px; color:' + colortrans+'">&nbsp&nbsp transmission : ' + transmission + '</label></b><br><br>';

		}
	}
	
	document.getElementById('idtitre').innerHTML = titre;
	
	var trio;
	
	var tabPosR;
	var tabRatio;
	var tabRX;
	var tabRY;
	
	var tabPosR_dad;
	var tabRatio_dad;
	var tabRX_dad;
	var tabRY_dad;
	
	var tabPosR_mom;
	var tabRatio_mom;
	var tabRX_mom;
	var tabRY_mom;
	
	var tabPosBA;
	var tabMother;
	var tabFather;
	var tabPatient;
	
	var tabplotRX=[];
   	var tabcnvRX=[];
   	var tabplotRY=[];
   	var tabcnvRY=[];
   	
   	var tabplotRX_dad=[];
   	var tabcnvRX_dad=[];
   	var tabplotRY_dad=[];
   	var tabcnvRY_dad=[];
   	
   	var tabplotRX_mom=[];
   	var tabcnvRX_mom=[];
   	var tabplotRY_mom=[];
   	var tabcnvRY_mom=[];
   	
   	var tabplotBAX=[];
   	var tabcnvBAX=[];
   	
   	var tabplotYM=[];
   	var tabcnvYM=[];
   	
   	var tabplotYF=[];
   	var tabcnvYF=[];
   	
   var tabplotYP=[];
   var tabcnvYP=[];
   
  	var graphDiv = document.getElementById('plot1');	
  	
	dataStorePlotValue.fetch({
			onComplete: function(items,request){
                	var item = items[0];
   
   					trio = parseInt(dataStorePlotValue.getValue(item, "PLOT"));
   			
   					 tabPosR = dataStorePlotValue.getValue(item, "POSratio"); 
   					 tabRatio = dataStorePlotValue.getValue(item, "RATIO");
  
   					tabRX = tabPosR.split(" ");
   					tabRY = tabRatio.split(" ");
 	
 					tabplotRX.push(0);
  					tabplotRY.push(0);
  					
   					for ( var i = 0 ; i < tabRX.length ; i++ ) {
   					
 						//points aberrants
   						if ( tabRY[i]  < -2 ) 
   						{
 								tabRY[i] = -2;
 						}
 						
 						if ( tabRY[i]  > 2) 
   						{
 								tabRY[i] = 2;
 						}
 						
  						if (  (parseInt(tabRX[i]) < parseInt(debcnv) ) || ( parseInt(tabRX[i]) > parseInt(fincnv) ) ) 
  						{
  							tabplotRX.push(tabRX[i]);
  							tabplotRY.push(tabRY[i]);
  						}
  						else
  						{
  							tabcnvRX.push(tabRX[i]);
  							tabcnvRY.push(tabRY[i]);
  						}
  					}
  		
  					var trace1 = {
  							x:tabplotRX,
  							y:tabplotRY,
  							yaxis: 'y1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
					};
						
					var trace2 = {
  							x:tabcnvRX,
  							y:tabcnvRY,
  							yaxis: 'y1',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'orange',
      								size:4
    						}
					};
					
				// ratio des parents
					
				if ( (trio == 1) || (trio ==3) )
				{	
						// 1 pour le pere
						tabPosR_dad = dataStorePlotValue.getValue(item, "POSratio_dad"); 
   					 	tabRatio_dad = dataStorePlotValue.getValue(item, "RATIO_dad");
  	
   						tabRX_dad = tabPosR_dad.split(" ");
   						tabRY_dad= tabRatio_dad.split(" ");
   						
   						tabplotRX_dad.push(0);
  						tabplotRY_dad.push(0);
  						tabplotRX_mom.push(0);
  						tabplotRY_mom.push(0);
 						
   						for ( var i = 0 ; i < tabRX_dad.length ; i++ ) {
   						
   							//points aberrants
   							if ( tabRY_dad[i]  < -2 ) 
   							{
 								tabRY_dad[i] = -2;
 							}
 						
 							if ( tabRY_dad[i]  > 2) 
   							{
 								tabRY_dad[i] = 2;
 							}
 
  							if (  (parseInt(tabRX_dad[i]) < parseInt(debcnv) ) || ( parseInt(tabRX_dad[i]) > parseInt(fincnv) ) ) 
  							{
  								tabplotRX_dad.push(tabRX_dad[i]);
  								tabplotRY_dad.push(tabRY_dad[i]);
  							}
  							else
  							{
  								tabcnvRX_dad.push(tabRX_dad[i]);
  								tabcnvRY_dad.push(tabRY_dad[i]);
  							}
  						}
  		
  						var trace3 = {
  							x:tabplotRX_dad,
  							y:tabplotRY_dad,
  							yaxis: 'y2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
						};
						
						var trace4 = {
  							x:tabcnvRX_dad,
  							y:tabcnvRY_dad,
  							yaxis: 'y2',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'lightblue',
      								size:4
    						}
						};
				}
					
				if 	( (trio == 2) || ( trio ==3) )
				{
						
						// pour la mere
						tabPosR_mom = dataStorePlotValue.getValue(item, "POSratio_mom"); 
   					 	tabRatio_mom = dataStorePlotValue.getValue(item, "RATIO_mom");
  
   						tabRX_mom = tabPosR_mom.split(" ");
   						tabRY_mom= tabRatio_mom.split(" ");
 	
   						for ( var i = 0 ; i < tabRX_mom.length ; i++ ) {
   						
   							//points aberrants
   							if ( tabRY_mom[i]  < -2 ) 
   							{
 								tabRY_mom[i] = -2;
 							}
 						
 							if ( tabRY_mom[i]  > 2) 
   							{
 								tabRY_mom[i] = 2;
 							}
 
  							if (  (parseInt(tabRX_mom[i]) <= parseInt(debcnv) ) || ( parseInt(tabRX_mom[i]) >= parseInt(fincnv) ) ) 
  							{
  								tabplotRX_mom.push(tabRX_mom[i]);
  								tabplotRY_mom.push(tabRY_mom[i]);
  							}
  							else
  							{
  								tabcnvRX_mom.push(tabRX_mom[i]);
  								tabcnvRY_mom.push(tabRY_mom[i]);
  							}
  						}
  		
  						var trace5 = {
  							x:tabplotRX_mom,
  							y:tabplotRY_mom,
  							yaxis: 'y3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'grey',
      								size:3
    						}
						};
						
						var trace6 = {
  							x:tabcnvRX_mom,
  							y:tabcnvRY_mom,
  							yaxis: 'y3',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'pink',
      								size:4
    						}
						};
						
				}
					
				if (trio == 3)
				{	
						tabPosBA = dataStorePlotValue.getValue(item, "POSball"); 
   						tabMother = dataStorePlotValue.getValue(item, "MOTHER");
 						tabFather = dataStorePlotValue.getValue(item, "FATHER");
 						
 				/*		MeanF = dataStorePlotValue.getValue(item, "MEANF");
 						MeanM = dataStorePlotValue.getValue(item, "MEANM");
 						
 						var label = 'mono allelic duplication';
 						if ( (MeanF < 50) &&  (MeanM <50) )
 						{
 						 	label = 'bi allelic duplication';
 						}
 				*/
   						var tabBAX = tabPosBA.split(" ");
   						var tabYM = tabMother.split(" ");
   						var tabYF = tabFather.split(" ");
   					
   						for ( var i = 0 ; i < tabBAX.length ; i++ ) {
   					
  							if (  (parseInt(tabBAX[i]) < parseInt(debcnv) ) || ( parseInt(tabBAX[i]) > parseInt(fincnv) ) ) 
  							{
  									tabplotBAX.push(tabBAX[i]);
  									tabplotYM.push(tabYM[i]);			
  									tabplotYF.push(tabYF[i]);
  							}
  							else
  							{
  									tabcnvBAX.push(tabBAX[i]);
  									tabcnvYM.push(tabYM[i]);
  									tabcnvYF.push(tabYF[i]);
  							}
  						}	
  	
  						var trace7 = {
  							x:tabplotBAX,
  							y:tabplotYM,
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'pink',
      								size:4
    						}
						};
				
						var trace8 = {
  							x:tabplotBAX,
  							y:tabplotYF,
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'lightblue',
      								size:4
    						}
						};
			
						var trace9 = {
  							x:tabcnvBAX,
  							y:tabcnvYM,
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'deeppink',
      								size:4
    						}
						};
				
						var trace10 = {
  							x:tabcnvBAX,
  							y:tabcnvYF,
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'blue',
      								size:4
    						}
						};
					
					
					
		/*				
						var trace11 = {
  							x:[0],
  							y:[MeanF],
  							yaxis: 'y4',
  							mode: 'markers',	
  							showlegend:false,
  							name:' ',
  							text: [label],
  							textfont : {
    								family:'Times New Roman'
  							},
  									
  							marker: {
      								color: 'blue',
      								symbol:['arrow-right-open'],
      								size:11,
    						},				
						};		
						
						var trace12 = {
  							x:[0],
  							y:[MeanM],
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'deeppink',
      								symbol:['arrow-right-open'],
      								size:11
    						},
						};
						
		     			if ( type == 1) // DEL
						{
							data=[trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9,trace10];
						}
						else
						{
							data=[trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9,trace10,trace11,trace12];
						}
		*/	
		
		
						data=[trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9,trace10];
						
						layout = {
  							yaxis1:   {domain: [0.0, 0.2], title : 'Wc ratio : Child'},
  							yaxis2:   {domain: [0.25, 0.45], title : ' father',color:'blue'},
  							yaxis3:   {domain: [0.5, 0.7], title : ' mother',color:'deeppink'},
  							yaxis4:	  {domain: [0.8, 1], title : 'Allelic balance', tickmode: "array", tickvals: [0, 25, 50, 75, 100], ticktext: ['0', '25', '50', '75', '100']},
  							xaxis : {title : 'genomic position'},
						};
			}
			else
			{
						// pour la balance allelique
						tabPosBA = dataStorePlotValue.getValue(item, "POSball"); 
   						tabPatient = dataStorePlotValue.getValue(item, "PATIENT");
 						
   						var tabBAX = tabPosBA.split(" ");
   						var tabYP = tabPatient.split(" ");
   						
   						for ( var i = 0 ; i < tabBAX.length ; i++ ) 
   						{
  							if (  (parseInt(tabBAX[i]) < parseInt(debcnv) ) || ( parseInt(tabBAX[i]) > parseInt(fincnv) ) ) 
  							{		
  									tabplotBAX.push(tabBAX[i]);
  									tabplotYP.push(tabYP[i]);								
  							}
  							else
  							{		
  									tabcnvBAX.push(tabBAX[i]);
  									tabcnvYP.push(tabYP[i]);
  							}
  						}	
  					
  						var trace7 = {
  							x:tabplotBAX,
  							y:tabplotYP,
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'lightgrey',
      								size:3
    						}
						};
						
						var trace8 = {
  							x:tabcnvBAX,
  							y:tabcnvYP,
  							yaxis: 'y4',
  							mode: 'markers',
  							showlegend:false,
  							name:' ',
  							marker: {
      								color: 'purple',
      								size:4
    						}
						};
						
						if (trio == 1) // on a que le pere
						{
								
								data=[trace1,trace2,trace3,trace4,trace7,trace8,];
								var layout = {
  									yaxis:   {domain: [0.1, 0.3], title : 'Wisecondor ratio'},
  									yaxis2:   {domain: [0.35, 0.55], title : ' father'},
  									yaxis4: {domain: [0.6, 1], title : 'Allelic frequencies',tickmode: "array", tickvals: [0, 25, 50, 75, 100], ticktext: ['0', '25', '50', '75', '100']},
  									xaxis : {title : 'Genomic positions'},
								};
						}
						
						if (trio == 2) // on a que la mere
						{
								
								data=[trace1,trace2,trace5,trace6,trace7,trace8,];
								var layout = {
  									yaxis:   {domain: [0.1, 0.3], title : 'Wisecondor ratio'},
  									yaxis3:   {domain: [0.35, 0.55], title : ' mother'},
  									yaxis4: {domain: [0.6, 1], title : 'Allelic frequencies',tickmode: "array", tickvals: [0, 25, 50, 75, 100], ticktext: ['0', '25', '50', '75', '100']},
  									xaxis : {title : 'Genomic positions'},
								};
						}
						
						if (trio == 0) // on a que l'enfant
						{	
							
							data=[trace1,trace2,trace7,trace8];
					
							var layout = {
  								yaxis:   {domain: [0.1, 0.5], title : 'Wisecondor ratio'},
  								yaxis4: {domain: [0.6, 1], title : 'Allelic frequencies',tickmode: "array", tickvals: [0, 25, 50, 75, 100], ticktext: ['0', '25', '50', '75', '100']},
  								xaxis : {title : 'Genomic positions'},
							};
						}
			}	
			
			var config = { modeBarButtonsToRemove: ['pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d','toggleSpikelines', 'hoverClosestCartesian', 'hoverCompareCartesian'], displayModeBar: true};
			Plotly.newPlot('plot1', data, layout,config);
			
			
			graphDiv.on('plotly_selected', function(eventData)
 				{
 
					var max = 0;
					var min = 250000000;
					
 					eventData.points.forEach(function(pt) {
  						min = Math.min(min,pt.x);
  						max = Math.max(max,pt.x);
 					});
 
 					var locus = chr+":"+min+"-"+max;
 					viewGenes(locus,project,patient);
 					
 					Plotly.restyle(graphDiv, data, layout,config);
				});
		}
  	});		
  	
}










