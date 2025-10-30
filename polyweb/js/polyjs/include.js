/**
 * @author pnitschk
 */

var url_path = "/cgi-bin/polymorphism-cgi/";


   function onSearchAll (input,name){
                var string = "'*"+input+"*'";
                var grid = dijit.byId(name);
                grid.queryOptions = {ignoreCase: true};
               grid.filter({complexQuery: "name:"+string+"OR description:"+string+"OR users:"+string+"OR patient_name:"+string+"OR capture_name:"+string+"OR new_pathogenic:"+string});
        };
        
        function onSearchAll_polyviewer (input,name){
            var string = "'*"+input+"*'";
            var grid = dijit.byId(name);
            grid.queryOptions = {ignoreCase: true};
           grid.filter({complexQuery: "name:"+string+"OR description:"+string+"OR users:"+string+" OR patient_name:"+string+" OR capture_name:"+string+"OR sites:"+string});
    };   
function param( name )
{
  name = name.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
  var regexS = "[\\?&]"+name+"=([^&#]*)";
  var regex = new RegExp( regexS );
  var results = regex.exec( window.location.href );

  if( results == null )
    return "";
  else
    return results[1];
}

function findColumnIndex(layout,name) {
		for (var i=0;i<layout.length;i++){
		
			if (layout[i]["name"] == name) return i; 
		}	
		alert("error can't find "+name +" for layout");
}

function findColumnIndex2(layout,name) {
		for (var i=0;i<layout.length;i++){
			for (var j = 0; j < layout[i].length; j++) {
				if (layout[i][j]["name"] == name) 
					return j;
			}
		}	
		alert("error can't find "+name +" for layout");
}	
//c'est ici que sont stockes toutes les url vers les cgi



var cache_ext="";
var href = window.location.href;

var url_track_ucsc = "http%3A%2F%2Fmendel.necker.fr%2Fcgi-bin%2Fpolymorphism-cgi%2Fucsc_bed.pl%3Fproject%3D";
if (href.indexOf("aptana", 1)>-1 ) {
	if (href.indexOf("cecile", 1) > -1) {
		url_path = "/cgi-bin/cfourrag/polymorphism-cgi/";
		url_track_ucsc = "http%3A%2F%2Fmendel.necker.fr%2Fcgi-bin%2Fpoly-dev%2Fpolymorphism-cgi%2Fucsc_bed.pl%3Fproject%3D";
	}
	else {
		//	cache_ext = "_dev";
		url_path = "/cgi-bin/pnitschk/polymorphism-cgi/";	
		url_track_ucsc = "http%3A%2F%2Fmendel.necker.fr%2Fcgi-bin%2Fpolymorphism-cgi%2Fucsc_bed.pl%3Fproject%3D";
	}
}

else if (href.indexOf("mbras", 1) > -1) {
	url_path = "/cgi-bin/mbras/polymorphism-cgi/";
}

else if (href.indexOf("archive", 1) > -1) {
	url_path = "/cgi-bin/archive/polymorphism-cgi/";
}

else if (href.indexOf("masson", 1) > -1) {
	url_path = "/cgi-bin/masson/polymorphism-cgi/";
}

else if (href.indexOf("eollivie", 1) > -1) {
	url_path = "/cgi-bin/eollivie/polymorphism-cgi/";
}
else if (href.indexOf("polydev", 1) > -1) {
	url_path = "/cgi-bin/polydev/polymorphism-cgi/";
}
else if (href.indexOf("preprod", 1) > -1) {
	url_path = "/cgi-bin/preprod/polymorphism-cgi/";
}
else if (href.indexOf("dev", 1) > -1) {
	url_path = "/cgi-bin/dev/polymorphism-cgi/";
}
// db or not
//var json_out = "/json_output";
var json_out = "/json_output_nodb";
var url_listProject = url_path+json_out+"/projects_json.pl";
var url_listProject_vector = url_path+json_out+"/projects_json_vector.pl";


//Filter url 

var filter_directory = "filter";

var url_view_filter = url_path+"/"+filter_directory+"/view_detail_filter.pl";
var url_view_filter_vector = url_path+"/"+filter_directory+"/view_detail_filter_vector.pl";
var url_region = url_path+"/"+filter_directory+"/load_regions.pl";
var url_delete_filter =  url_path+"/"+filter_directory+"/delete_filter.pl";
var url_get_filter =  url_path+"/"+filter_directory+"/get_filter.pl";
var url_save_filter =  url_path+"/"+filter_directory+"/save_filter.pl";



var url_report_variations =  url_path+"/xls_output/variations_report.pl";

var url_auth = url_path+"authentification.pl";
var url_query =  url_path+"projectList_json.pl";
var url_insert =  url_path+"insertData.pl";
var url_chart = url_path+json_out+"/chart_json.pl";
var url_annot =  url_path+json_out+"/viewer_json.pl";
var url_cover =  url_path+"viewCover.pl";
var url_cache = url_path+"cache_json.pl";
var url_igv = url_path+json_out+"/load_igv.pl";
var url_find_name =  url_path+json_out+"/search_name.pl";
var url_region_genes = url_path+json_out+"/regions_genes_json.pl";
var url_cache_variation = url_path+json_out+"/variations_json.pl";
var url_cache_variation_vector = url_path+json_out+"/variations_json_vector.pl";
var url_cache_variation_excel_vector = url_path+json_out+"/polyquery.pl";
var url_cache_transcripts = url_path+json_out+"/transcripts_json.pl";
var url_cache_patient = url_path+json_out+"/patients_json.pl";
var url_stat_polyweb =url_path+json_out+"/stats_polyweb.pl";

var is_cache =  url_path+"isCache.pl";
//var url_cache_graph = "/cgi-bin/cache_graph.pl";
var url_report = url_path+"report_nomenclature.pl";
var url_layout = url_path+json_out+"/send_layout.pl";
var url_layout_vector = url_path+json_out+"/send_layout_vector.pl";
var url_cache_gene =  url_path+json_out+"/genes_json.pl";
var url_cache_gene_vector =  url_path + 'json_output_nodb/polyquery.pl';

var url_valid = url_path+"/validation_variation/validations.pl";
var url_valid_acmg = url_path+"/validation_variation/validations_acmg.pl";
var url_track_ucsc = escape("http://mendel.necker.fr/"+url_path+json_out+"/ucsc_bed.pl?project=");
var url_genes_infos = url_path+json_out+"/genes_infos.pl";
var url_genes_infos_vector = url_path+json_out+"/interface_gene_json.pl";

var url_follow = url_path+json_out+"/follow_projects.pl";

// validation and coverage project url
var url_check_project_update = url_path+"/validation_variation/check_project_update.pl";
var url_patients = url_path+"/validation_variation/patients_json.pl";
var url_validation_genes = url_path+"/validation_variation/genes_json.pl";
var url_coverage = url_path+"/validation_variation/coverage_exons.pl";
var url_primer_coverage = url_path+"/validation_variation/coverage_primers.pl";
var url_save_exons = url_path+"/validation_variation/exons_todo.pl";
var url_todo =  url_path+"/validation_variation/todo_json.pl";
var url_report_patient =  url_path+"/validation_variation/patient_report.pl";

var url_report_patient_cache =  url_path+"/validation_variation/variations_editor.pl";
var url_update_validations =  url_path+"/validation_variation/update_validations.pl";
var url_chart_transcript = url_path+"/validation_variation/transcript_chart.pl";
var url_chart_exon = url_path+"/validation_variation/exon_chart.pl";
var url_chart_exon_hc = url_path+"/validation_variation/exon_chart_hc.pl";
var url_chart_gene = url_path+"/validation_variation/gene_chart.pl";
var url_chart_primer= url_path+"/validation_variation/primer_chart.pl";
var url_summary = url_path+"/validation_variation/summary.pl";
var url_summary_panel = url_path+"/validation_variation/summary_panel.pl";
var url_header = url_path+"/validation_variation/header.pl";
var url_diag_tree = url_path+"/validation_variation/bundle_json.pl";
var url_list_transcripts = url_path+"/validation_variation/list_transcripts_json.pl";	
var url_coverage_json = url_path+"/validation_variation/coverage_json.pl";  
//var url_coverage_overview = url_path+"/validation_variation/table_image.pl";
var url_coverage_overview = url_path+"/validation_variation/table_image.PROD.pl";
var url_coverage_overview_defidiag = url_path+"/validation_variation/table_image_updated.pl";
var url_primer_coverage_overview = url_path+"/validation_variation/table_image_primer_coverage.pl"; 
var url_cnv_overview = url_path+"/validation_variation/table_image_dude.pl"; 
//var url_cnv_overview = url_path+"/validation_variation/table_images_cnv.PROD.pl";
var url_cached_list_primers = url_path+"/validation_variation/preload_list_primer.pl"; 
var url_cached_coverage = url_path+"/validation_variation/preload_coverage.pl"; 
var url_cnv = url_path+"/validation_variation/table_dude.pl";
//var url_cnv = url_path+"/validation_variation/image_cnv.PROD.pl";
var url_problem = url_path+"/problem_capture/report_problem.pl";
var url_polydiag = url_path+"/validation_variation/patient_report.pl";
var url_polyviewer = url_path+"/validation_variation/patient_report.pl";
var url_hgmd_view = url_path+"/hgmd/view.pl";
var url_report_quality = url_path+"/report/quality.pl";
var url_igv_js = url_path+"/igv/load_igv.pl";
var url_CNV_Patient = url_path + "/manta/Retrieve_allSV_Patient.pl";
// mes tests
// var url_parserManta = url_path + "/manta/parser_VcfManta.pl";
