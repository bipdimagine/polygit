require(["dojox/charting/widget/Legend", "dojox/charting/action2d/Tooltip","dojox/charting/Chart","dojox/charting/DataChart",
 "dojox/charting/axis2d/Default", "dojox/charting/plot2d/Markers","dojox/charting/plot2d/Columns","dojox/charting/plot2d/ClusteredColumns", "dojox/charting/themes/PurpleRain" , "dojo/ready"]);
 
 /* Chart display */
var type_chart = 1;
 function load_graph_transcript (pn,tr,ex){
 	  dijit.byId("graph_exon").show();
           hideShowWait();
        dojo.byId("chart_exon").innerHTML ="";
//     dojo.byId("waiting_chart").innerHTML ="<img src='/icons/Polyicons/wait_big.gif'>";
    
   // dojo.byId("title").innerHTML = "<center>"+pn+" "+tr+" "+ex+"</center>";
    var url = url_chart_transcript+"?project="+projectName+"&patients="+pn+"&transcripts="+tr;
    
        var xhrArgs = {
                          url: url,
                          handleAs: "json",
                          preventCache: true,
                          obj:this,
                         load: function(data){
                          // Replace newlines with nice HTML tags.
                            hideShowWait();
                         dijit.byId("graph_exon").show();
                              show_graph_bar(data.items,'Coverage exons/patient  :'+pn);
                 },
                 error: function(error){
                 	   hideShowWait();
                     this.innerHTML = "An unexpected error occurred: " + error;
                   }   
            };
          
                    // Call the asynchronous xhrGet
                   var deferred = dojo.xhrGet(xhrArgs);
    // type_chart= 2;
    // var store = new dojo.data.ItemFileWriteStore({
                // url:url
        // });
        // store.fetch({
                // onComplete:show_graph_exon_ori
            // });
}
 
 

 
function load_graph_exon (pn,tr,ex,start,end,raw){
	   dijit.byId("graph_exon").show();
	    hideShowWait();
	    if (dijit.byId('waiting')) {
                    dijit.byId('waiting').show();
                }   
    
      //   dijit.byId("popup_coverage").close();
         //dojo.byId("waiting_chart").innerHTML ="<img src='/icons/Polyicons/wait_big.gif'>";
        dojo.byId("chart_exon").innerHTML ="";
//        dojo.byId("title").innerHTML = "<center>"+tr+" "+ex+"</center>";
        
          var padding = 20;
     /*   if (padding_spinner){
        padding = padding_spinner.value;
        }
        */
   var url = url_chart_exon_hc+"?project="+projectName+"&patients="+pn+"&transcripts="+tr+"&exon="+ex+"&span="+padding+"&start="+start+"&end="+end+"&type="+raw;
        var xhrArgs = {
                          url: url,
                          handleAs: "json",
                          preventCache: true,
                          obj:this,
                         load: function(data){
                          // Replace newlines with nice HTML tags.
                      	    hideShowWait();
                              show_graph_exon(data.items);
                 },
                 error: function(error){
                     this.innerHTML = "An unexpected error occurred: " + error;
                   }   
            };
          
                    // Call the asynchronous xhrGet
                         var deferred = dojo.xhrGet(xhrArgs);

   
}
var nb_bar;
var col_names = [];
function show_graph_bar (items,title){
	var tt ="";
	tt = tt+items[0].title;
	if (tt.length >0){
		 dojo.byId("title_graph").innerHTML =tt;
	}
	else {
		dojo.byId("title_graph").innerHTML =title;
	}
    
  col_names = [];
  var colorSets = ['#F45B5B', '#90ED7D', '#8AE234'];
       
  
   var cat = new Array();
   var mean = new Array();
      var mean_all = new Array();
      var points = new Array();
      nb_bar =0;
      	points[0] = ["0","0","0"];
      	var csv;
      	csv ="x,coverage,project\n";
     for (i=0;i<items.length;i++){
     	cat[i] = items[i].seq;
     	mean[i] = items[i].cov;
     	mean_all[i] = items[i].covm;
     	col_names[i] = items[i].seq;
     	csv +=items[i].seq+","+items[i].cov+","+items[i].covm+"\n";
     	points[i] = [i,items[i].cov,items[i].covm];
     }
     
//      dojo.byId("waiting_chart").innerHTML ="";
    dojo.byId("chart_exon").innerHTML ="";
g = new Dygraph(
              document.getElementById("chart_exon"),
              points,
              {
               labelsDiv: document.getElementById('status'),
              	labels:["Patients","coverage","mean"],
                legend: 'always',
              //  title: 'Daily Widget Sales',
                includeZero: true,
               axes: {
                x: {
                  axisLabelFormatter: function(d, gran, opts) {
                  	 if (col_names[d] == undefined ) return "" ; 
                      return col_names[d];
                  },
                  pixelsPerLabel: 20,
                  axisLabelFontSize :6,
                },
               
               },
                 highlightCallback: function(e, x, pts, row) {
                	pts_infos2(e,x,pts,items);
              }, 
                colors: colorSets,
                animatedZooms: false,
                //drawXGrid: false,
                plotter: barChartPlotter,
                dateWindow: [-1, items.length],
                drawXGrid:false,
              }
          ); 
}

function pts_infos2(e,x,pts,items){
	var s = document.getElementById("status");
	
		
	
	var i = pts[0].xval;
	var seq = items[i].seq;
	var mean =  items[i].cov;
	var mean_project = items[i].covm;
	    s.innerHTML = "<b>"+seq+": "+mean+"  Project: "+mean_project+"</b>";
}


function darkenColor(colorStr) {
        // Defined in dygraph-utils.js
        var color = Dygraph.toRGB_(colorStr);
        color.r = Math.floor((255 + color.r) / 2);
        color.g = Math.floor((255 + color.g) / 2);
        color.b = Math.floor((255 + color.b) / 2);
        return 'rgb(' + color.r + ',' + color.g + ',' + color.b + ')';
      }

function barChartPlotter(e) {
        var ctx = e.drawingContext;
        var points = e.points;
        var y_bottom = e.dygraph.toDomYCoord(0);

       
		nb_bar ++;
        // Find the minimum separation between x-values.
        // This determines the bar width.
        var min_sep = Infinity;
        for (var i = 1; i < points.length; i++) {
          var sep = points[i].canvasx - points[i - 1].canvasx;
          if (sep < min_sep) min_sep = sep;
        }
        var size = 5;
         var	color  = "rgba(169, 255, 150,1.0)";
        if (nb_bar % 2 == 0 ){
        	size = 10;
        	color  = "rgba(244, 91, 91,1.0)";
        }
        var bar_width = Math.floor(2.0 / size * min_sep);
	
        // Do the actual plotting.
         var size = 2;
        for (var i = 0; i < points.length; i++) {
          var p = points[i];
         
          var center_x = p.canvasx;
		ctx.fillStyle  = "rgba(200, 200, 200, 1.0)";
			 ctx.fillRect(center_x - bar_width / size +1, p.canvasy +1,
              bar_width+1, y_bottom - p.canvasy);
              ctx.fillStyle = e.color; 
          ctx.fillRect(center_x - bar_width / size, p.canvasy,
              bar_width, y_bottom - p.canvasy);
  				
  	
          ctx.strokeRect(center_x - bar_width / 2, p.canvasy,
              bar_width, y_bottom - p.canvasy);
              
        }
      }
      

function show_graph_bar2 (items,title){
   var t = items.data;
   var cat = new Array();
   var mean = new Array();
      var mean_all = new Array();
      
     for (i=0;i<items.length;i++){
     	cat[i] = items[i].seq;
     	mean[i] = items[i].cov;
     	mean_all[i] = items[i].covm;
     }
      dojo.byId("waiting_chart").innerHTML ="";
    dojo.byId("chart_exon").innerHTML ="";
     var chart = new Highcharts.Chart({
     chart: {
            type: 'column',
             renderTo: 'chart_exon',
         
        },
        title: {
            text: title
        },
        subtitle: {
            text: ''
        },
        colors: ['#f45b5b', '#90ed7d', '#90ed7d', '#f7a35c', '#8085e9', 
   '#f15c80', '#e4d354', '#2b908f', '#f45b5b', '#91e8e1'],
        xAxis: {
            categories: cat,
            crosshair: true
        },
        yAxis: {
            min: 0,
            title: {
                text: ''
            }
        },
        tooltip: {
            headerFormat: '<span style="font-size:10px">{point.key}toto</span><table>',
            pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
                '<td style="padding:0"><b>{point.y:.1f} </b></td></tr>',
            footerFormat: '</table>',
            shared: true,
            useHTML: true
        },
        plotOptions: {
            column: {
                grouping: false,
                shadow: true,
                borderWidth: 0
            }
        },
        series: [ 
        {
            name: 'project',
            data: mean_all,
              pointPadding: 0.3,
            pointPlacement: -0.1
           

        },
        {
            name: 'coverage',
            data: mean,
           pointPadding: 0.4,
           pointPlacement: -0.1
			
        },
        
       ]
        
    });
 
}

 
function load_graph_one_exon_all_patients (tr,ex,start,end){
	
        dijit.byId("graph_exon").show();
           hideShowWait();
//         dojo.byId("waiting_chart").innerHTML ="<img src='/icons/Polyicons/wait_big.gif'>";
        dojo.byId("chart_exon").innerHTML ="";

    var url = url_chart_gene+"?project="+projectName+"&transcripts="+tr+"&exon="+ex+"&start="+start+"&end"+end;
        type_chart= 2;
       var xhrArgs = {
                          url: url,
                          handleAs: "json",
                          preventCache: true,
                          obj:this,
                         load: function(data){
                          // Replace newlines with nice HTML tags.
                         		hideShowWait();
                              show_graph_bar(data.items,'coverage patients/exon : '+tr+':'+ex);
                 },
                 error: function(error){
                     this.innerHTML = "An unexpected error occurred: " + error;
                   }   
            };
          
                    // Call the asynchronous xhrGet
                         var deferred = dojo.xhrGet(xhrArgs);

    
}
function show_graph_exon (items){

   //dojo.byId("waiting_chart").innerHTML ="";
    dojo.byId("chart_exon").innerHTML ="";
     dijit.byId("graph_exon").show();
      dojo.byId("title_graph").innerHTML=items.title;
   
    var exonic = items.exonic;
    var utr =  items.utr;
    var intronic =  items.intronic;
  
   g = new Dygraph(
    document.getElementById("chart_exon"),
    items.csv, // path to CSV file
    {
    		//animatedZooms: true,
    	  	//title:items.title,
    	  	strokeWidth: 2,
    	  	showRangeSelector: true,
    	  	labels: items.labels,
    	  	 ylabel: 'Coverage',
    	  	  legend: 'always',
    	  	  axes: {
    	  	  	Y: {
    	  	  		max:3000,	
    	  	  	},
    	  	  },
    	  	 'Coverage': {
                  pointSize: 1.5,
                },
             /*   'Mean Coverage': {
                  pointSize: 1.5
                },*/
         highlightCallback: function(e, x, pts, row) {
                	pts_infos(e,x,pts,items);
              
              },       
         underlayCallback: function(canvas, area, g) {
         	for (var i=0;i<exonic.length;i++){
              var bottom_left = g.toDomCoords(exonic[i].from, -20);
              var top_right = g.toDomCoords(exonic[i].to, +20);

              var left = bottom_left[0];
              var right = top_right[0];
              canvas.fillStyle = "rgba(223, 240, 216, 1.0)";
              canvas.fillRect(left, area.y, right - left, area.h);
             }
             for (var i=0;i<intronic.length;i++){
              var bottom_left = g.toDomCoords(intronic[i].from, -20);
              var top_right = g.toDomCoords(intronic[i].to, +20);

              var left = bottom_left[0];
              var right = top_right[0];
              canvas.fillStyle = "rgba(230, 226, 235,1.0)";
              canvas.fillRect(left, area.y, right - left, area.h);
             }
             for (var i=0;i<utr.length;i++){
              var bottom_left = g.toDomCoords(utr[i].from, -20);
              var top_right = g.toDomCoords(utr[i].to, +20);

              var left = bottom_left[0];
              var right = top_right[0];
              canvas.fillStyle = "rgba(250, 250, 250,1.0)";
              canvas.fillRect(left, area.y, right - left, area.h);
             }
            },
            
          
             
            
            }          // options
  );
 
}



function pts_infos(e,x,pts,items){
	var s = document.getElementById("status");
		var seq = items.sequences[x];
		var pose =  x-items.exon_start;
	var genomic = ((items.genomic_start*1)+(x*1));

	if (pose == 0){
		pose =1;
	}
	
	else if (x<items.exon_start){
		pose = x-items.exon_start ;
	}
	else if (x>items.exon_end){
		pose = "+"+(x-items.exon_end );
	}

	var exon_position = x - items.splice5;
	if (exon_position > items.exon_length ){
		exon_position = "+"+(exon_position-items.exon_length);
	} 
	    s.innerHTML = x + "<b>  Genomic:&nbsp"+genomic+"</b> -- exon Position : " +pose+ "--  sequences :"+seq+"";
}

function show_graph_exon2 (items){
   var t = items.data;
	console.log(items);
        dojo.byId("waiting_chart").innerHTML ="";
    dojo.byId("chart_exon").innerHTML ="";
     var chart = new Highcharts.Chart({
    chart: {
            renderTo: 'chart_exon',
            zoomType: 'xy'
        },
        title: {
           // text: items.gene_name+"-"+items.transcript_name+"-"+items.exon_name,
            text: items.title,
        },
        legend: {
            layout: 'vertical',
            align: 'right',
            verticalAlign: 'middle',
            borderWidth: 1
        },
        subtitle: {
            text: document.ontouchstart === undefined ?
                    'Click and drag in the plot area to zoom in' :
                    'Pinch the chart to zoom in'
        },
        xAxis: {
           plotBands: items.legende,
            type: 'numeric',
            minRange: 50 // fourteen days,
           
        },
        yAxis: {
            title: {
                text: 'Depth',
                type:'numeric',
                /* min: 10 ,
                 max:1500,*/
            },
            min:0,
        },
  
        plotOptions: {
            area: {
                fillColor: {
                    linearGradient: { x1: 0, y1: 0, x2: 0, y2: 1},
                    stops: [
                        [0, Highcharts.getOptions().colors[2]],
                        [1, Highcharts.Color(Highcharts.getOptions().colors[2]).setOpacity(0).get('rgba')]
                    ]
                },
                marker: {
                    radius: 2
                },
                lineWidth: 1,
                states: {
                    hover: {
                        lineWidth: 1
                    }
                },
                threshold: null,
                  turboThreshold: 10000,
            }
        },

        series: [
        {
            type: 'line',
            name: items.names[0],
            pointInterval: 1,
             pointStart: 0,
           // pointStart: items.start*-1,
            data: t[0]["data"],
        },
         {
            type: 'line',
            name: items.names[1],
            pointInterval: 1,
                 pointStart: 0,
           // pointStart:items.start*-1,
            data: t[1]["data"],
        }]
        
    });
 
}

 
function load_graph_gene (tr){
        dijit.byId("graph_exon").show();
           hideShowWait();
        dojo.byId("chart_exon").innerHTML ="";
        var pn="";
    var url = url_chart_gene+"?project="+projectName+"&patients="+pn+"&transcripts="+tr;
        type_chart= 2;
    // var store = new dojo.data.ItemFileWriteStore({
                // url:url
        // });
        // store.fetch({
                // onComplete:show_graph_bar
            // });

        var xhrArgs = {
                          url: url,
                          handleAs: "json",
                          preventCache: true,
                          obj:this,
                         load: function(data){
                          // Replace newlines with nice HTML tags.
                         hideShowWait();
                              show_graph_bar(data.items,'coverage patients/transcript :'+tr);
                 },
                 error: function(error){
                 	   hideShowWait();
                     this.innerHTML = "An unexpected error occurred: " + error;
                   }   
            };
          
                    // Call the asynchronous xhrGet
                         var deferred = dojo.xhrGet(xhrArgs);
    
}




function show_graph_exon_ori (item){
    var t = new Array();
    var z = new Array();
    var m = new Array; 
    for (i=0;i<item.length;i++){
    
        var x = new Array();
            var s = item[i]["seq"];
            t[i] = item[i]["cov"]*1.0;
            m[i] = item[i]["covm"]*1.0;
        //alert(s);
        z[i] = {value:i+1, text:s+""};
    }
    
     dojo.byId("waiting_chart").innerHTML ="";
    dojo.byId("chart_exon").innerHTML ="";
  var c = new dojox.charting.Chart("chart_exon");
  if (type_chart == 2){
      c.addPlot("default", {type: dojox.charting.plot2d.ClusteredColumns,gap:2});
     }
     else {
         c.addPlot("default", {type: dojox.charting.plot2d.Markers,gap:2});
     }
      c.addAxis("x", {
    labels: z, dropLabels: false
    });
       
    c.addAxis("y", {vertical: true, fixLower: "major", fixUpper: "major"});
        c.setTheme(dojox.charting.themes.PurpleRain);
        c.addSeries("patient", t);
        c.addSeries("Project Mean", m);
      new dojox.charting.action2d.Tooltip(c, "default", {
       text: function(o){
          return z[o.index].text+" <br>coverage :  "+t[o.index]+" <br> project mean: "+m[o.index]+"";
       }
    });
    
        c.render();
   
 var legend = new dojox.charting.widget.Legend({chart: c, horizontal: true}, "legend");
    
}