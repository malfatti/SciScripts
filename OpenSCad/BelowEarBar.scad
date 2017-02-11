BarWid = 5; BarLen = 45; BarFit = 5; CircleRes = 50;

module RoundHolder() {
    translate([0, 0, BarWid/2]) rotate([90, 0, 0]) difference() {
        cylinder(r=BarWid/2, h=BarLen, $fn=CircleRes, center=true);
        cylinder(r=BarWid/2-1, h=BarLen, $fn=CircleRes, center=true);
        translate([0, 0, BarLen/2-2]) cylinder(r=BarWid/2-0.5, h=BarFit, $fn=CircleRes, center=true);
        translate([0, 0, -BarLen/2+2]) cylinder(r=BarWid/2-0.5, h=BarFit, $fn=CircleRes, center=true);
    }
}

module SqHolder() {
    translate([0, 0, BarWid/2]) rotate([90, 0, 0]) difference() {
        cube([BarWid, BarWid, BarLen], true);
        cube([BarWid-2, BarWid-2, BarLen], true);
        translate([0, 0, BarLen/2-2]) cube([BarWid-1, BarWid-1, BarFit], true);
        translate([0, 0, -BarLen/2+2]) cube([BarWid-1, BarWid-1, BarFit], true);
    }
}

module ArcCircle() {
    difference() {
        cylinder(r=BarWid*2+2, h=2, $fn=CircleRes);
        cylinder(r=BarWid*2, h=2, $fn=CircleRes);
    }
}

module QuarterArc() {
    intersection() {
       translate([BarWid*2+2, BarWid*2+2, 0]) cube([(BarWid*2+2)*2, (BarWid*2+2)*2, 2], true);
        ArcCircle();
    }
}

module HalfArc() {
    intersection() {
        translate([0, BarWid*2+2, 0]) cube([(BarWid*2+2)*2, (BarWid*2+2)*2, 2], true);
        ArcCircle();
    }
}

module RoundJoint() {
    union() {
        translate([0, BarWid*2+0.5, BarWid/2]) rotate([90, 0, 0]) difference() {
            cylinder(r=BarWid/2, h=3, $fn=CircleRes, center=true);
            cylinder(r=BarWid/2-1, h=3, $fn=CircleRes, center=true);
        }
        translate([0, BarWid*2+4.5, BarWid/2]) rotate([90, 0, 0]) difference() {
            cylinder(r=BarWid/2-0.5, h=BarFit, $fn=CircleRes, center=true);
            cylinder(r=BarWid/2-1, h=BarFit, $fn=CircleRes, center=true);
        }
    }
}

module SqJoint() {
    union() {
        translate([0, BarWid*2+0.5, BarWid/2]) rotate([90, 0, 0]) difference() {
            cube([BarWid, BarWid, 3], true);
            cube([BarWid-2, BarWid-2, 3], true);
        }
        translate([0, BarWid*2+4.5, BarWid/2]) rotate([90, 0, 0]) difference() {
            cube([BarWid-1, BarWid-1, BarFit], true);
            cube([BarWid-2, BarWid-2, BarFit], true);
        }
    }
}

module ExtJoint() {
    union() {
        translate([0, BarWid*2+0.5, BarWid/2]) rotate([90, 0, 0]) difference() {
            cylinder(r=BarWid/2, h=3, $fn=CircleRes, center=true);
            translate([0, 0, 3]) cylinder(r=BarWid/2-1, h=6, $fn=CircleRes, center=true);
        }
        translate([0, BarWid*2+4, BarWid/2]) rotate([90, 0, 0]) difference() {
            cylinder(r=BarWid+2, h=BarFit, $fn=CircleRes, center=true);
            translate([0, 0, -1]) cylinder(r=BarWid, h=BarFit-1, $fn=CircleRes, center=true);
            translate([0, 0, 3]) cylinder(r=BarWid/2-1, h=6, $fn=CircleRes, center=true);
        }
    }
}

module Arc() {
    union() {
        
        HalfArc();
        rotate([90, 0, 0]) HalfArc();
        difference() {
            rotate([0, -90, 0]) QuarterArc();
            translate([0, BarWid*2+2, BarWid/2]) rotate([90, 0, 0]) cylinder(r=BarWid/2-1, h=8, $fn=CircleRes, center=true);
        }
    }
}

module Tip() {
    union() {
        translate([0, 1.5, 0]) rotate([90, 0, 0]) difference() {
            cylinder(r=BarWid/2, h=3, $fn=CircleRes, center=true);
            cylinder(r=BarWid/2-1, h=3, $fn=CircleRes, center=true);
        }
        translate([0, 7, 0]) rotate([90, 0, 0]) difference() {
            cylinder(r2=BarWid, r1=2, h=BarFit+3, $fn=CircleRes, center=true);
//            translate([0, 0, -1]) cylinder(r=BarWid, h=BarFit-1, $fn=CircleRes, center=true);
            translate([0, 0, 0]) cylinder(r=BarWid/2-1, h=8, $fn=CircleRes, center=true);
        }
    }
}

module RoundBar() {
    union() {
        rotate([0, 0, 180]) translate([0, -20, 0]) Arc();
        rotate([0, 0, 180]) translate([-0.5, -20, 0]) RoundJoint();
        translate([0, -BarLen/2, 0]) RoundHolder();
    }
}

module SqBar() {
    union() {
        rotate([0, 0, 180]) translate([0, -20, 0]) Arc();
        rotate([0, 0, 180]) translate([-0.5, -20, 0]) SqJoint();
        translate([0, -BarLen/2, 0]) SqHolder();
    }
}

module ExtBar() {
    union() {
        rotate([0, 0, 180]) translate([0, -BarWid*2-2-4.5, 0]) Arc();
        rotate([0, 0, 180]) translate([-0.5, -BarWid*2-2-4.5, 0]) ExtJoint();
    }
}

translate([15, 50, 0]) SqBar();
rotate([180, 180, 0]) translate([-35, -30, 0]) SqBar();
translate([55, 50, 0]) RoundBar();
rotate([180, 180, 0]) translate([-75, -30, 0]) RoundBar();
rotate([90, 0, 0]) translate([95, 0, -50]) ExtBar();
rotate([90, 0, 0]) translate([95, 0, -70]) ExtBar();
rotate([270, 0, 0]) translate([85, -11, 30]) Tip();
rotate([270, 0, 0]) translate([105, -11, 30]) Tip();