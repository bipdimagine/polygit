/**
 * @author pnitschk
 */

dojo.declare("transcript", null, {
		//acts like a java constructor
	
		constructor: function(x,y,w,name,color){
		this.w = w;	
		this.x = x;	
		this.y= y;
		this.name = name;
		this.color = color;
		this.exon = new Array();
		this.h_text=8;
	},
	
	
	drawText: function (){	
		
	/*	makeText(surface_transText, {
                x: 0,
                y:this.y+space*1.7+heightrectRef,
                text: this.name,
                align: "start"
            }, {
                family: "Verdana",
                size: "7pt",
                weight: "bold"
            }, "black");
            */
		
	},
	drawGraph:function(){
		this.group = surface_trans.createGroup();
		drawRect(surface_trans,1,this.y,dataRef.length*echx-1,sizeRectTranscript+4+this.h_text,this.color,decx,1);
		drawRect(surface_trans,this.x*echx,this.y+3,this.w*echx,sizeRectTranscript-4,[11,11,12],decx,1);
		 this.g = surface_trans.createGroup();
		 var start_text = this.exon[0].x*echx;
		//this.g.connect("onmouseover", dojo.hitch(this, "mousein"));
		
		 
		
		
		 for (i in this.exon) {
		 	var ex = this.exon[i];
			  	//var color1 = [0,255,80,1];
				var color1 = [96,255,70,1];
				var color2 = [128,128,128,1];
        
            	if (ex.strand == -1){
                	color1 = [255,27,99,1];
					color2 = [128,128,128,1];
            	}
			drawRect(this.g,ex.x*echx,this.y,ex.w*echx,sizeRectTranscript+2,color1,decx);
			
			if (ex.utrstart){
					
					drawRect(this.g,ex.utrstart*echx,this.y,(ex.utrend-ex.utrstart)*echx,sizeRectTranscript+2,color2,decx);
			}
		 }
		 
		 makeText(surface_trans, {
                x: start_text,
                y:this.y+space*1.7+this.h_text,
                text: this.name,
                align: "start"
            }, {
                family: "Verdana",
                size: "6pt",
                weight: "bold"
            }, "black");
	},
	addExon:function(x,w,strand,utrstart,utrend){
	
		this.exon.push({
			x: x,
			w: w,
			strand: strand,
			utrstart: utrstart,
			utrend: utrend
		});
		
	},
	mousein : function (e) {
			document.getElementById("over").innerHTML= "Genes : "+this.name;
	},
	
	});