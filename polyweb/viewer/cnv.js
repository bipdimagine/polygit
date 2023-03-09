/**
 * @author pnitschk
 */


dojo.declare("cnv", [variation], {
		//acts like a java constructor
	
	draw: function() 
	{
		this.ishide = false;
		this.group = surface.createGroup();
		this.grouphighlight = surface.createGroup();
		var xstart = this.x*echx;
		var xend = this.data.length*echx;
		for ( i in this.patient){
			var color1 = [0,0,0,0.01];
		var ystart = this.patient[i].y;
			drawRect(this.group,xstart,ystart,xend,6,color1,decx);
		 
		}	
			this.group.connect("onmouseover", dojo.hitch(this, "mousein"));
	 		this.group.connect("onmouseout", dojo.hitch(this, "mouseout"));
	
			
	},
	highlight: function(){
		
		
		var text="";
		for ( i in this.patient){	
			text += this.patient[i].name + " ";		
		 	groupPatients[this.patient[i].name].textHiglight();
		}
		
		document.getElementById("over").innerHTML="Patient: "+text;
		
		
	},
	redraw :function (){
		document.getElementById("over").innerHTML="";
		for (i in this.patient) {
			groupPatients[this.patient[i].name].writeName();
		}
		
	},
});	  