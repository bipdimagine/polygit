
// TODO: step 1 - import IGV dans page HTML   ->  <script src="https://cdn.jsdelivr.net/npm/igv@2.11.2/dist/igv.min.js"></script>
// TODO: step 2 - import igv_view_methods.js  ->  <script type="text/javascript" src="../js/polyjs/igv_view_methods.js" djConfig="parseOnLoad: true"></script>
// TODO: step 3 - creation Dialog/Div         ->  <div dojoType="dijit.Dialog" jsId="dialog_igv" id="dialog_igv" title="IGV View" style="background-color:white;overflow: auto;" > <div class="container-fluid" id="div_igv" style="height:90%;width:90%;"></div> </div>
// TODO: step 4 - appel methode               ->  view_web_igv_bam("dialog_igv", "div_igv", locus, bams, names, genome)


var url_gene_bed, url_fasta, url_cyto;
var is_waiting = false;

function view_web_igv_bam_simple(div_name, locus, file, name, genome) {
	var array_bam, array_name;
	if   (file.match(';')) { array_bam = file.split(";"); }
	else                   { array_bam = file.split(","); }
	if   (name.match(';')) { array_name = name.split(";"); }
	else                   { array_name = name.split(","); }
	
	var div = document.getElementById(div_name);
	div.innerHTML = "";

	if (genome){
		if (url_gene_bed) {}
		else { url_gene_bed = "/public-data/"+genome+"/igv/gencode.gtf.gz"; }
		if (url_fasta) {}
		else { url_fasta = "/public-data/"+genome+"/genome/fasta/all.fa"; }
		if (url_cyto) {}
		else { url_cyto = "/public-data/"+genome+"/igv/cytoband.txt"; }
	} 
	else {
		genome = 'HG19';
		const regex = /HG38/g;
		if (file.match(regex)) {
			genome = 'HG38';
		}
		if (url_gene_bed) {}
		else { url_gene_bed = "/public-data/"+genome+"/igv/gencode.gtf.gz"; }
		if (url_fasta) {}
		else { url_fasta = "/public-data/"+genome+"/genome/fasta/all.fa"; }
		if (url_cyto) {}
		else { url_cyto = "/public-data/"+genome+"/igv/cytoband.txt"; }
	}
	
	var list_bams = [];
	var list_tracks = [];
	var track = {};
	track['name'] = "Genes";
    track['type'] = "annotation";
    track['format'] = "gtf";
    track['sourceType'] = "file";
    track['url'] = url_gene_bed;
    track['indexURL'] = url_gene_bed+".tbi";
    track['order'] = Number.MAX_VALUE;
    track['visibilityWindow'] = 300000000;
    track['displayMode'] = "EXPANDED";
    track['showSoftClips'] = true;
	list_tracks.push(track);
	for (var i=0;i<array_bam.length;i++){
		var track = {};
		track['type'] = 'alignment';
		track['format'] = 'bam';
		track['name'] = array_name[i] + ' BAM';
		track['url'] = window.location.origin + array_bam[i];
		track['indexURL'] = window.location.origin + array_bam[i] + '.bai';
		track['autoHeight'] = true;
    	track['samplingDepth'] = 100.;
        track['coverageThreshold'] = 1;
        track['coverageQualityWeight'] = true;
        track['visibilityWindow'] = 1000000;
        track['colorBy'] = "pairOrientation";
        track['showSoftClips'] = true;
		list_tracks.push(track);
		list_bams.push(window.location.origin + array_bam[i]);
	}

    if (typeof locus === 'undefined') { locus = 'chr1:1-249,250,621'; }
    var options = {
        showNavigation: true,
        showRuler: true,
        showAllChromosomes:true,
        reference: {
            fastaURL: url_fasta,
            cytobandURL: url_cyto,
        },
        locus: locus,
        trackDefaults: {
            bam: {
                coverageThreshold: 0.2,
                coverageQualityWeight: true,
                height: 300,
                maxHeight: 10000,
                visibilityWindow: 1000000,
                showSoftClips: true,
                
                
            }
        },
        palette: [
            ["#00A0B0", "#6A4A3C", "#CC333F", "#EB6841"]
        ],
        tracks: list_tracks,
        trackDefaults: {
            bam: {
                coverageThreshold: 0.2,
                coverageQualityWeight: true,
                height: 3000,
                maxHeight: 10000,
                visibilityWindow: 1000000,
                showSoftClips: true,
            },
            bed: {
                height: 30,
                maxHeight: 100,
                visibilityWindow: 1000000,
            }
        }
    };
    
    setTimeout(function() {
		igv.createBrowser(div, options).then(function (browser) {
		    browser.search(locus);
		    var tracks = list_bams.join(',');
		    if (url_fasta.match(/HG19/)) {
		    	tracks = tracks + ',' + window.location.origin + '/' + url_gene_bed;
		    	launch_igv_tool('', tracks, locus);
		    }
		    else if (url_fasta) {
		    	url_fasta = window.location.origin + '/' + url_fasta;
		    	launch_igv_tool(url_fasta, tracks, locus);
		    }
		    else { launch_igv_tool('', tracks, locus); }
		    is_waiting = false;
		});
    }, 500);
}

var browser;
var browser_igv;

function view_web_igv_bam(dialog_name, div_name, locus, file, name) {
	const regex = /HG38/g;
	var genome = 'hg19';
	if (file.match(regex)) {
		genome = 'hg38';
	}
	
	file = file.replaceAll(",", ";");
	name = name.replaceAll(",", ";");
	
	if (is_waiting == true) { return; } 
	is_waiting = true;
	
	var array_bam, array_name;
	if   (file.match(';')) { array_bam = file.split(";"); }
	else                   { array_bam = file.split(","); }
	if   (name.match(';')) { array_name = name.split(";"); }
	else                   { array_name = name.split(","); }
	
	var div = document.getElementById(div_name);
	div.innerHTML = "";
	if (dialog_name == '') {}
	else {
		resise_div_igv(div_name);
		setTimeout(function() { dijit.byId(dialog_name).show(); }, 200);
	}
	
    if (typeof locus === 'undefined') { locus = 'chr1:1-249,250,621'; }
    
    var list_tracks = [];
    var list_names = name.split(';');
    var list_files = file.split(';');
	for (i=0; i<list_files.length; i++) {
		var this_file = list_files[i];
		
		var index = this_file + '.bai';
	    var type_file = 'bam';
	    if (this_file.match(".cram")) {
	    	type_file = 'cram';
	    	index = this_file + '.crai';
	    }
	    
	    var track = {};
	    track['type'] = 'alignment';
	    track['format'] = type_file;
	    track['url'] = this_file;
	    track['indexURL'] = index;
	    track['name'] = list_names[i];
	    track['height'] = 600;
        track['showSoftClips'] = true;
	    list_tracks.push(track);
	}
    
    var options =
        {
            genome: genome,
            locus: locus,
            tracks:list_tracks
        }
    
    setTimeout(function() {
    	if (browser_igv == 1) {
    		document.getElementById("div_igv").remove();
    		div = document.createElement('div');
    		div.setAttribute("id", "div_igv");
			div.classList.add('container-fluid');
			div.style.margin = '20px';
			document.getElementById(dialog_name).appendChild(div);
    	}
   		browser = igv.createBrowser(div, options);
   		browser_igv = 1;
		is_waiting = false;
    }, 500);
}

/*
function getTrackGenesAnnotations() {
	var track = {};
	track['type'] = 'annotation';
	track['format'] = 'bed';
	track['name'] = 'HG19 Genes';
	track['url'] = "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz";
	track['indexURL'] = "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/genes/refGene.hg19.bed.gz.tbi";
	track['order'] = Number.MAX_VALUE;
	track['visibilityWindow'] = 300000000;
    track['displayMode'] = "EXPANDED";
    track['showSoftClips'] = true;
    return track;
}

function getTrackPhase3WgsVariants() {
	var track = {};
	track['type'] = 'variant';
	track['format'] = 'vcf';
	track['name'] = 'Phase 3 WGS Variants';
    track['showSoftClips'] = true;
	track['url'] = "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz";
	track['indexURL'] = "https://s3.amazonaws.com/1000genomes/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi";
    return track;
}
*/
function resise_div_igv(div_name) {
    var body = document.body;
    var html = document.documentElement;
    var height = Math.max( body.scrollHeight, body.offsetHeight, html.clientHeight, html.scrollHeight, html.offsetHeight ) - 50;
    var width  = Math.max( body.scrollWidth, body.offsetWidth, html.clientWidth, html.scrollWidth, html.offsetWidth ) - 150;
    $("#" + div_name).height(height);
    $("#" + div_name).width(width);
}


function igv_go_to(locus) {
	is_waiting = true;
    igv.browser.search(locus);
    is_waiting = false;
}
            
