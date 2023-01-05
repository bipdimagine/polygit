/**
 * @author pnitschk
 */


dojo.declare("variation", null, {
		//acts like a java constructor
		constructor: function(x,y,color,data){
		
		this.x =x;
		this.y = y-2;
		this.color = color;
		this.data = data;
		this.patient = new Array();
		this.w = data.length;
	},
	draw: function() 
	{
		this.ishide = false;
		this.group = surface_pat.createGroup();
		this.grouphighlight = surface_pat.createGroup();
		
		for ( i in this.patient){
 	var aShape = this.group.createImage({width:8, height: 8,src:"/icons/Polyicons/led"+this.color+".png"}).
 	applyTransform({dx: this.x*echx+decx, dy: this.patient[i].y}) ;
	
		 aShape.connect("onmouseover", dojo.hitch(this, "mousein"));
		 
		 
		}	
	
		if (this.data.validate == 1){
			 var aShape2 = this.group.createImage({width: 10, height: 10,src:"/icons/Polyicons/kgpg.png"}).
			applyTransform({dx:  this.x*echx+decx, dy: legend_height});
		}
		if (this.data.validate == -1){
			 var aShape2 = this.group.createImage({width: 10, height: 10,src:"/icons/Polyicons/kgpg.png"}).
			applyTransform({dx:  this.x*echx+decx, dy: legend_height});
		}
		
			
	},
	setColor: function(color){
		this.color = color;
		
	},
	addPatient: function(name,y){
	
		this.patient.push({y:y,name:name});
		
	},
	highlight: function(){
		
		surface_pat.remove(this.group);
		if (this.group2) {
			surfaceMap.remove(this.group2);
		}
		this.group = surface_pat.createGroup();
		this.group2 = surfaceMap.createGroup();
		var text= "";
		for ( i in this.patient){
		text += this.patient[i].name + " ";
		 groupPatients[this.patient[i].name].textHiglight();
		  var aShape = this.group.createImage({width: 6, height: 6,src:"/icons/Polyicons/ledyellow.png"}).
		applyTransform({dx: this.x*echx+decx, dy: this.patient[i].y}) ;
		 aShape.connect("onclick", dojo.hitch(this, "mouseclick"));
		 aShape.connect("onmouseout", dojo.hitch(this, "mouseout"));
		}
		
		document.getElementById("over").innerHTML= "Patients: "+text;
	
		//this.group.clear();
		
	},
	
	size1 : function(){
		return (this.sizeIndicator);
	},
	redraw :function (){
		document.getElementById("over").innerHTML="";
		for (i in this.patient) {
			
			groupPatients[this.patient[i].name].writeName();
		}
		surface_pat.remove(this.group);
		surfaceMap.remove(this.group2);
		this.draw();
	},
	hide : function (){
		this.ishide = true;
		surface_pat.remove(this.group);
	},
	translate : function (dx){
		this.group.applyTransform({dx: -1*dx});
	},
	mouseclick : function(){
		//	var store = readStore(varurl);
	/*		store.fetch({
				query: {
					name: this.data.name
				},
				onComplete: detailVar
			
			});
		*/	
		//detailVar(projectName,this.data.name,this.data.reference,this.data.contig);
//		var fenetre =parent.window.open("detailVar.html",name,config=' toolbar=no, menubar=no, scrollbars=yes, resizable=no, location=no, directories=no, status=no');		
		this.data.view =1;
		this.redraw();
			
	},
	mousein : function (e) {
		if (lasthighlight){
			lasthighlight.redraw();
			lasthighlight= null;
		}
		this.highlight();
		lasthighlight = this;
	},
	mouseout : function () {
		this.redraw();
		lasthighlight= null;
		
	},
	isUTR : function (){
		if (this.data.type == "UTR") return true; 
		return false;
	},
	isCoding : function () {
		if (this.data.type == "Coding") return true; 
		return false;
	},
	isEnsembl : function () {
		if (this.data.ensembl == 1) return true; 
		return false;
	},
	isIntronic : function () {
		if (this.data.type == "Intronic" || this.data.type == "Intergenic") return true;  
		return false;
	},
	hideorshow : function (){
	
		if (this.data.filterscore == "weak" && !showWeak){
			this.hide();
			}
		else if (this.data.filterscore == "medium" && !showMedium){
			this.hide();
		}
		else if (this.data.filterscore == "strong" && !showStrong){
			this.hide();
		}
		else if (this.isEnsembl() && !showEnsembl){
			this.hide();
			
		}
		else if (this.isCoding() && !showCoding){
			this.hide();
		}
		else if (this.isUTR() && !showUTR){
			this.hide();
		}
		else if (this.isIntronic() && !showIntronic){
			this.hide();
		}
		else {
			if (this.data.filterscore == "weak") statistics2["nbweak"]++;
			else if (this.data.filterscore == "medium") statistics2["nbmedium"]++;
			else if (this.data.filterscore == "strong") statistics2["nbstrong"]++;
			if (this.isCoding())  statistics2["nbcoding"]++;
			else if (this.isIntronic())  statistics2["nbintronic"]++;
			else if (this.isUTR())  statistics2["nbutr"]++;
			if (this.isEnsembl())  statistics2["nbdbsnp"]++;
			if (this.data.validate == 1 )  statistics2["nbvalidate"]++;
			 if (this.ishide) this.redraw();
			 
		} 
		
	}
	
	
});	  
	 