BarWid = 5; BarLen = 45; CircleRes = 500;

module RoundHolder() {
    translate([0, BarWid/2, 0]) cylinder(r=BarWid/2, h=BarLen, $fn=CircleRes, center=true);
}

module SqHolder() {
    translate([0, BarWid/2, 0]) cube([BarWid, BarWid, BarLen], true);
}

module ArcCircle() {
    difference() {
        translate([0, BarWid, 0]) cylinder(r=BarWid*2+2, h=2, $fn=CircleRes);
        translate([0, BarWid, 0]) cylinder(r=BarWid*2, h=2, $fn=CircleRes);
    }
}

module VertArc() {
    intersection() {
       rotate([0, 90, 0]) ArcCircle();
       translate([1, BarWid*2+2, -BarWid*2]) rotate([0, 90, 0]) cube([(BarWid*2+2)*2, (BarWid*2+2)*2, 2], true); 
    }
}
module FrontArc() {
    intersection() {
        translate([-BarWid*2-1, 0, 0]) cube([(BarWid*2+2)*2, (BarWid*2+2)*2, 2]);
        translate([1, 0, 0]) ArcCircle();
    }
}
module Arc() {
    VertArc();
    FrontArc();
}
//translate([-15-1, 0, BarLen/2+BarWid*2]) Arc();
//translate([-15, 0, 0]) RoundHolder(); 
//translate([-5, 0, 0]) RoundHolder(); 
//translate([5, 0, 0]) SqHolder();
//translate([15, 0, 0]) SqHolder();
translate([0, 0, BarWid*2+2]) Arc();

/*
module Holder(){
    module Outer(){
        union(){
            cube([HolderSide+Thick, HolderSide+Thick, Thick+Thick]);
            translate([OutCupCenter, OutCupCenter, Thick]) cylinder(r=OutCupDiam/2, h=CupLen, $fn=CircleRes);
        }
    }

    module Inner(){
        union(){
            translate([HolderCenter, HolderCenter, 0]) cube([HolderSide, HolderSide, Thick]);
            translate([OutCupCenter, OutCupCenter, Thick]) cylinder(r=CupDiam/2, h=CupLen, $fn=CircleRes);
        }
    }
    
    difference(){Outer(); Inner();}
}

module Holes(){
    translate([Hole1Center, Hole1Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    translate([Hole2Center, Hole2Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    translate([Hole1Center, Hole2Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    translate([Hole2Center, Hole1Center, Thick]) cube([HoleSide, HoleSide, Thick]);
    
    translate([Screw1X, Thick/2, Screw1Z]) rotate([90,0,0]) cylinder(r=ScrewDiam/2, h=Thick/2, $fn=CircleRes);
    translate([Screw1X, Screw2Y, Screw2Z]) rotate([90,0,0]) cylinder(r=ScrewDiam/2, h=Thick, $fn=CircleRes);
    
}
difference(){Holder(); Holes();}
*/
