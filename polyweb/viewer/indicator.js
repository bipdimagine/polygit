/**
 * @author pnitschk
 */
	/* Object Indicator */

dojo.declare("Indicator", null, {
		//acts like a java constructor
		constructor: function(surface,echxMap,lref){
		this.surface = surface;
		this.size = (lref*echxMap);
		this.height = 20;
		this.click = false;
		this.lastx = decx;
		
	},
	draw: function(zoom) 
	{
		this.zoom = zoom;	
		this.sizeIndicator = this.size/zoom;
		this.group = this.surface.createGroup();
	
		if (zoom==1) return ;
		this.group.connect("onmousedown",  dojo.hitch(this,"handleMouseDown"));		
		this.group.connect( "onmouseup",  dojo.hitch(this,"handleMouseUp"));
	
		this.group.createRect({
			x: decx,
			y: 2,
			height: this.height-1,
			width: this.sizeIndicator-3
		})
		.setFill([0,0,255,0.4])
		.setStroke({ color: "grey", width: 1});	
		
		this.group.createRect({
			x: decx,
			y: 2,
			height: this.height-3,
			width: this.sizeIndicator-3
		})	
		.setStroke({ color: "black", width: 2});	
	},
	
	size1 : function(){
		return (this.sizeIndicator);
	},
	
	redraw :function (zoom){
		this.surface.remove(this.group);
		this.draw(zoom);		
	},
	
	translate : function (dx){
		this.group.applyTransform({dx: -1*dx});
	},
	
	mousein : function(){
		this.inside=true;
	},
	
	move : function (x){
		var dx = (this.lastx - x);
		this.lastx= x;
		if (previousx  <dx){
			dx =previousx;			
		}
		translate(dx);
		previousx -= dx;
	},
	
	handleMouseDown : function (event){
        container_position = container_transcripts_position;  
		this.lastx = event.clientX - container_position.x,
		this.click=true;
		dojo.stopEvent(event);	
	},
	
	handleMouseMove: function (event){
		if(!this.click) return;
        container_position = container_transcripts_position;  
		var x = event.clientX - container_position.x;	
		this.move(x);
		dojo.stopEvent(event);
	},
	
	handleMouseUp: function (event){
		dojo.stopEvent(event);
	}
});
