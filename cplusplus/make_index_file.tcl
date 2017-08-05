set traj [lindex $argv 0]
set psf [lindex $argv 1]
mol new $psf type psf
mol addfile $traj type xtc waitfor all

proc write_file {selection} {
	set sel [atomselect top "${selection}"]
	set ind [$sel get index]

	set MIN [lindex [lsort -integer -increasing $ind] 0]
	set MAX [lindex [lsort -integer -decreasing $ind] 0]

	set wat [atomselect top "water"]
	set w_ind [$wat get index]
	set first_w [lindex [lsort -integer -increasing $w_ind] 0]
	set last_w [lindex [lsort -integer -decreasing $w_ind] 0]

	set nump [$sel num]
	set noh [atomselect top "$selection and noh"]
	set nump_noh [$noh num]

	puts "(selection) (first protein index) (last protein index) (number of protein atoms) (number of protein heavy atoms) (first water) (last water)"
	puts "$selection $MIN $MAX $nump $nump_noh $first_w $last_w"

	$sel delete
	set fil [open "indices.txt" w]
	for {set i 0} {$i<[llength $ind]} {incr i} {
		set I [lindex $ind $i]
		set sel_I [atomselect top "index $I"]
		puts $fil [format {%8s%8d} "[string range [$sel_I get type] 0 0]" "$I"]
		$sel_I delete
	}
	close $fil
}

write_file "resname TAS"
#write_file "protein"

exit
