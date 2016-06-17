ScrewD = 6; ScrewRes = 100;
Height = 42; Thick = 5; SqBase = 15;

module Base(){
    union(){
        cube([17.5, SqBase, Thick]);
        cube([Thick, SqBase, Height]);
        translate([0, 0, 25]) cube([69.5, SqBase, Thick]);
        translate([52, 0, 0]) cube([17.5, SqBase, Thick]);
        translate([66.5, 0, 0]) cube([Thick, SqBase, Height]);
    }
}

module Holes(){
    translate([10, 7.5, 0]) cylinder(h=Thick, r=ScrewD/2, $fn=ScrewRes);
    translate([59.5, 7.5, 0]) cylinder(h=Thick, r=ScrewD/2, $fn=ScrewRes);
    translate([0, 7.5, 35]) rotate([0, 90, 0]) cylinder(h=Thick, r=ScrewD/2, $fn=ScrewRes);
    translate([66.5, 7.5, 35]) rotate([0, 90, 0]) cylinder(h=Thick, r=ScrewD/2, $fn=ScrewRes);
}

difference(){Base(); Holes();}