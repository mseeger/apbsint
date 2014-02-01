# Runs through all files *.h,*.cc,*.c in this directory and its subdirectories
# (one level only) and exchanges expression $source against $target.
# The files are overwritten. Make a copy first!

# User-def. variables:
set numExch 1
#set source(1) {#define [A-Z]+_H\n}
#set target(1) "\\0\n#if HAVE_CONFIG_H\n#  include <config.h>\n#endif\n"
#set source(1) {#define [A-Z]+_[A-Z]+_H\n}
#set target(1) "\\0\n#if HAVE_CONFIG_H\n#  include <config.h>\n#endif\n"
#set source(1) { \* Module: [A-Za-z]+\n}
#set target(1) " * Library source file\n\\0"
set source(1) {\&errstr}
set target(1) "errstr"

# Make file list
#set fileList [list]
set dirNames [glob -nocomplain */]
#set dirNames [glob -nocomplain lhotse/*/]
#append dirNames " [glob -nocomplain src/*/]"
#set dirNames [glob -nocomplain src/*/]
#set fileList [glob "lhotse/\{*.cc,*.h,*.c\}"]
#append fileList " [glob "src/\{*.cc,*.h,*.c\}"]"
#set fileList [glob "src/\{*.cc,*.h,*.c\}"]
set fileList [glob -nocomplain "\{*.cc,*.h,*.c\}"]
foreach dir $dirNames {
    append fileList " [glob -nocomplain "\{${dir}*.cc,${dir}*.h,${dir}*.c\}"]"
}

# Do substitutions: print names of all files which are changed
foreach fname $fileList {
    set fid [open $fname r]
    set text [read -nonewline $fid]
    close $fid
    set numCh 0
    for {set i 1} {$i<=$numExch} {incr i} {
	set numCh [expr $numCh + [regsub -all -- $source($i) $text $target($i) text]]
    }
    if {$numCh>0} {
	puts $fname
	set fid [open $fname w]
	puts $fid $text
	close $fid
    }
}
