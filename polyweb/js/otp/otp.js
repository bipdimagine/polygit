const inputs = document.getElementById("inputs_otp"); 
  
inputs.addEventListener("input_otp", function (e) { 
    const target = e.target; 
    const val = target.value; 
  	
    if (isNaN(val)) { 
        target.value = ""; 
        return; 
    } 
  
    if (val != "") { 
        const next = target.nextElementSibling; 
        if (next) { 
            next.focus(); 
        } 
    } 
}); 
  
inputs.addEventListener("keyup", function (e) { 
    const target = e.target; 
    const key = e.key.toLowerCase(); 
	const regex = /[a-zA-Z]/g;
	if (String(key).match(regex)) {
	  	if (key == "arrowup") { return; }
	  	if (key == "arrowdown") { return; }
	  	if (key == "arrowleft") {
	        const previous = target.previousElementSibling; 
	        if (previous) { 
	            previous.focus(); 
	        } 
	  		return;
	  	}
	  	if (key == "arrowright") {
	        const next = target.nextElementSibling; 
	        if (next) { 
	            next.focus(); 
	        } 
	  		return;
	  	}
	    if (key == "backspace" || key == "delete") { 
	        target.value = ""; 
	        const prev = target.previousElementSibling; 
	        if (prev) { 
	            prev.focus(); 
            	prev.value = ""; 
	        } 
	        return; 
	    }
        target.value = ""; 
	    return;
    }
    else {
        const next = target.nextElementSibling; 
        if (next) { 
            next.focus(); 
            next.value = ""; 
        } 
        return; 
    }
});

function okfunction_otp() {
	var my_otp_verif = document.getElementById('passwd_otp').value;
}

function checkPasswordOtp(username){
    dojo.style("otp_login", "visibility", "hidden");
	var my_code = document.getElementById('input_otp_1').value;
	my_code += document.getElementById('input_otp_2').value;
	my_code += document.getElementById('input_otp_3').value;
	my_code += document.getElementById('input_otp_4').value;
	my_code += document.getElementById('input_otp_5').value;
	my_code += document.getElementById('input_otp_6').value;
    if (!my_code) {
    	otp_auth();
    }
    else {
        var url1 = url_listProject+"?action=otp&mycode="+my_code+"&login="+username;
        $.getJSON( url1, function( data ) {
            var verif = data.otp_verif;
        	if (verif == 'yes') {
            	var d = new Date();
				d.setTime(d.getTime() + 12*60*60*1000);
				//d.setTime(d.getTime() + 60);
				var expires = "expires=" + d.toGMTString();
        		document.cookie = "otp=ok; SameSite=strict; Secure;" + expires;
        		dojo.style("otp_login", "visibility", "hidden");
        		okfunction_view_case();
        	}
        	else if (verif == 'no') {
        		location.reload();
        		exit(0);
        	}
        	else { location.reload(); }
        	return;
	    })
	}
}

var last_otp_code;
function otp_auth(otp_code) {
	if(otp_code) { last_otp_code = otp_code; }
	document.getElementById("span_otpcode").innerHTML = last_otp_code;
	dojo.style("otp_login", "visibility", "visible");
}

function checkCookieOtp() {
	var verif = dojo.cookie("otp");
	if (verif && verif == 'ok') {
		return 1;
	}
	return;
}


function onpaste_input_otp(e) {
	let paste = (event.clipboardData || window.clipboardData).getData("text");
	var listTmp = String(paste).split('');
	if (listTmp.length == 6) {
		document.getElementById('input_otp_1').value = listTmp[0];
		document.getElementById('input_otp_2').value = listTmp[1];
		document.getElementById('input_otp_3').value = listTmp[2];
		document.getElementById('input_otp_4').value = listTmp[3];
		document.getElementById('input_otp_5').value = listTmp[4];
		document.getElementById('input_otp_6').value = listTmp[5];
	}
}


