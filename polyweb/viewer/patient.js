/**
 * @author pnitschk
 */


dojo.declare("patient1", null, {
		//acts like a java constructor
	
		constructor: function(y,w,name,color){				
		
		
		this.w = w;
		this.y= y;
		this.color = color;
		this.name = name;
		this.trace = new Array();
		this.traceforward = new Array();
		this.rect = {};
	}
	,
change_coverage : function(){
	for (z in this.rect) {
		for (i in this.rect[z]) {
			
			this.rect[z][i].hide();
		}
	}
	var sel = "cov5";
	if (select_coverage) sel = select_coverage.get("value");
	for (i in this.rect[sel]) {
		
		this.rect[sel][i].draw();
	}
	
	
	
	},	
mouseclick : function(event){
	for (z in this.rect) {
		for (i in this.rect[z]) {
			this.rect[z][i].hide();
		}
	}
	/* var x = event.clientX - container_position.x -3;
	 var pos= parseInt((x/echx)+indicator.lastx/echxMap+genomic_start);
	 if (chr_number) {
	 	displayInIGV(chr_number.attr('value'), pos, pos);
	 }*/
	},
	drawText: function (){	
	
	/*this.groupText = surface_patText.createGroup();
		 drawRect(this.groupText,0,this.y+heightrectRef-1,100,sizeRect+3,this.color,0,1);
                    makeText(this.groupText, {
                        x: 0,
                        y: this.y + space * 1.5 + heightrectRef,
                        text: this.name,
                        align: "start"
                    }, {
                        family: "Verdana",
                        size: "6pt",
                    
                    }, "black");
		*/
	},
	
	drawGraph:function(){
		
		var rect =  drawRect(surface_pat,1,this.y,this.w*echx-1,sizeRect+3,this.color,decx,1);
		rect.connect("onclick", dojo.hitch(this, "mouseclick"));
		var color = {};
		color["cov0"] = [255, 0, 0, 0.9];
		color["cov5"] = [0, 120, 0, 0.9];
		color["cov10"] = [0,175, 0, 0.9];
		color["cov15"] = [0, 225, 0, 0.9];
		color["cov20"] = [0, 255, 0, 0.9];	
		 for (i in this.trace) {
		 	var tr = this.trace[i];
			var zz= tr.strand;
			
			  if (!this.rect[zz]) this.rect[zz] = new Array() ;
			  this.rect[zz][i] = new iRect(tr.x,this.y+1, tr.w,sizeRect,surface_pat,color[zz]);
			  this.rect[zz][i].draw() ;	
				
			
			
		 }
		 this.writeName();
		this.change_coverage();
		 // this.groupText2 = surface_pat.createGroup();
		
	},
	textHiglight:function(){
		surface_pat.remove(this.groupText2);
		 this.groupText2 = surface_pat.createGroup();
		 
		  this.addName("red");				
		
		},
	addTrace:function(x,y,strand){
		this.trace.push({
			x: x,
			w: (y-x),
			strand: strand
		});
		
	},
	writeName:function(){

	if (this.groupText2) 	surface_pat.remove(this.groupText2);
	this.groupText2 = surface_pat.createGroup();
	
	 this.addName("#000000");
	
		
	},
	addName:function(color){
		
		 makeText(this.groupText2, {
                        x: 10,
                        y: this.y+ space * 1.5 +5,
                        text: this.name,
                        align: "start"
                    }, {
                        family: "Verdana",
                        size: "7pt",
						weight: "bold"
                    },color);
    
				
		
								
			makeText(this.groupText2, {
                        x: this.w*echx-100,
                        y: this.y+ space * 1.5 +5,
                        text: this.name,
                        align: "start"
                    }, {
                        family: "Verdana",
                        size: "8pt",
						weight: "bold"
                    },color);
		
	}
	
});


dojo.declare("iRect", null, {
		//acts like a java constructor
		constructor: function(x,y,w,h,surface_pat,color){
		//(tr.x,this.y+sizeTrace+1, tr.w,surface_pat,[255, 0, 0, 0.8]);
		this.x =x;
		this.y = y;
		this.color = color;
		this.h = h;
		this.w = w;
		this.surface = surface_pat;
	},
	draw: function() 
	{
		
	this.rect =  drawRect(this.surface, this.x * echx, this.y, this.w * echx, this.h,this.color, decx,1);	
	this.rect.connect("onclick", dojo.hitch(this, "mouseclick"));
	
			
	},
	setColor: function(color){
		this.color = color;
		
	},
	addPatient: function(name,y){
	
		this.patient.push({y:y,name:name});
		
	},
	highlight: function(){
	
		
	},
	
	size1 : function(){
		
	},
	redraw :function (){
		
		this.draw();
	},
	hide : function (){
		this.ishide = true;
		surface_pat.remove(this.rect);
	},
	
	mouseclick : function(){
	
	 if (chr_number) {
	 	displayInIGV(chr_number.attr('value'), this.x+genomic_start,this.x+this.w+genomic_start);
	 }
	 else {
	 	displayInIGV(chromosome_name, this.x+genomic_start,this.x+this.w+genomic_start);
	 }
	
	},
	
	
	
});	  
