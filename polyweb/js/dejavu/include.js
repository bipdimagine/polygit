/**
 * @author mbras
 */

var href = window.location.href;

var url_path = "/cgi-bin/polymorphism-cgi/polydejavu";

if (href.indexOf("mbras", 1) > -1) {
	url_path = "/cgi-bin/mbras";
}
if (href.indexOf("aptana", 1) > -1) {
    url_path = "/cgi-bin/pnitschk/polymorphism-cgi/polydejavu";
}