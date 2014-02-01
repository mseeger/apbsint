# Runs through all files *.h,*.cc,*.c in this directory and its subdirectories
# (one level only) and greps for command line argument.

set dirNames [glob */]
set what [lindex $argv 0]

puts "DIRECTORY: / ---------------------------------------------------"
foreach fname [glob -nocomplain "\{*.cc,*.h,*.c\}"] {
    set flag [catch {exec grep $what $fname} outp]
    if {$flag==0} {
	puts "\[$fname\]: $outp"
    }
}
foreach dir $dirNames {
    puts "DIRECTORY: $dir ---------------------------------------------------"
    foreach fname [glob -nocomplain "\{${dir}*.cc,${dir}*.h,${dir}*.c\}"] {
	set flag [catch {exec grep $what $fname} outp]
	if {$flag==0} {
	    puts "\[$fname\]: $outp"
	}
    }
}
