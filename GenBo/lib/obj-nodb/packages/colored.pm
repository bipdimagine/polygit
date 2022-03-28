package colored;
use strict;
use Term::ANSIColor;

 my $stabilo_color = {red => ['black ON_BRIGHT_RED'],
 					 			  green => ['black ON_BRIGHT_GREEN'],
 					 			   yellow => ['black ON_BRIGHT_YELLOW'],
 					 			    cyan => ['black on_bright_cyan'],
 					 			    blue => ['black on_bright_blue'],
 					 			    magenta => ['black on_bright_magenta'],
 					 			    white => ['black on_bright_white'],
 };
 
 my $color_text = {red => ['red'],
 					 			  green => ['green'],
 					 			   yellow => ['yellow'],
 					 			    cyan => ['cyan'],
 					 			    blue => ['blue'],
 					 			    magenta => ['magenta'],
 					 			    white => ['white'],
 					 			    black =>['black'],
 };
 
 sub stabilo {
 	my ($color,$text,$string) =@_;
 	$color="white" unless exists  $stabilo_color->{$color};
 	if ($string ){
 		return text($stabilo_color->{$color},$text);
 	}
 	print1($stabilo_color->{$color},$text);
 	
 }
 
 sub print1 {
 	my ($color,$text) = @_;
 	 print  colored $color,$text  ;
		print  color 'reset';
		print  "\n";
 }
 sub text {
 	my ($color,$text) = @_;
 	return  colored $color, $text; 
 }
 sub print_color{
 	my ($color,$text,$string) =@_;
 	$color="black" unless exists  $color_text->{$color};
 	if ($string ){
 		return text($color_text->{$color},$text);
 	}
 	print1 ($color_text->{$color},$text);
 	
 }
 
 1;
 