"""

"""

from avlGeo import avlGeo
import xml.etree.ElementTree as ET
import math

######################
vsp_filepath = 'C:/Users/dalto/Documents/GrabCAD/Advanced 2021-2022/Analysis/OpenVSP/SB2_prelim_2.vsp3'
avl_out_path = 'SB2_prelim_2.avl'

control_surface_tolerance = 0.02  # does not create new control surface section if close enough to existing section
digit_precision = 3
######################

avl = avlGeo(avl_out_path, True)

######################
avl.mach = 0

avl.IYsym = 0
avl.IZsym = 0
avl.Zsym = 0

avl.xref = 1.05
avl.yref = 0
avl.zref = 0

Cspace = 'cosine'
Sspace = 'cosine'

Vortex_density_span = 5  # number of vortices to place per vsp dimension
Vortex_density_chord = 8

#Vortex_density_span = 10/12  # number of vortices to place per vsp dimension
#Vortex_density_chord = 10/12
######################

tree = ET.parse(vsp_filepath)
root = tree.getroot()


class ControlSurfaceSec:
    def __init__(self, name, u_loc, xhinge):
        self.name = name
        self.u_loc = u_loc
        self.xhinge = xhinge
        self.is_written = False


def define_surface(geom, totalSpan, chord):
    nchordwise = round(chord * Vortex_density_chord)
    nspanwise = round(totalSpan * Vortex_density_span)

    name = geom.find('ParmContainer/Name').text

    xloc = float(geom.find('ParmContainer/XForm/X_Location').get('Value'))
    yloc = float(geom.find('ParmContainer/XForm/Y_Location').get('Value'))
    yloc_rel = float(geom.find('ParmContainer/XForm/Y_Rel_Location').get('Value'))
    zloc = float(geom.find('ParmContainer/XForm/Z_Location').get('Value'))
    xrot = float(geom.find('ParmContainer/XForm/X_Rotation').get('Value'))
    
    if float(geom.find('ParmContainer/Sym/Sym_Planar_Flag').get('Value')) == 2:  # check symmetry and set yduplicate accordingly
        if float(geom.find('ParmContainer/Sym/Sym_Ancestor_Origin_Flag').get('Value')) == 0:  # check for "Object" button (ancestor origin)
            avl.addSurface(name, nchordwise, 0, Cspace, nspanwise, Sspace, round(yloc, digit_precision))
        else:
            avl.addSurface(name, nchordwise, 0, Cspace, nspanwise, Sspace, round(yloc-yloc_rel, digit_precision))
    else:
        avl.addSurface(name, nchordwise, 0, Cspace, nspanwise, Sspace)

    relative_twist = bool(float(geom.find('ParmContainer/WingGeom/RelativeTwistFlag').get('Value')))
    relative_dihedral = bool(float(geom.find('ParmContainer/WingGeom/RelativeDihedralFlag').get('Value')))

    control_surfaces = []

    num_sections = len(geom.findall('WingGeom/XSecSurf/XSec')) + 1

    for control_surface in geom.findall('Geom/SubSurfaces/SubSurface/SubSurfaceInfo/[Type=\'3\']/..'):
        # sets u locations as a percent of total (not projected) span
        u_start = (float(control_surface.find('ParmContainer/SS_Control/UStart').get('Value')) - 1 / num_sections) / (
                (num_sections - 2) / num_sections)  # weird transform needed due to how openvsp handles subsurfaces
        u_end = (float(control_surface.find('ParmContainer/SS_Control/UEnd').get('Value')) - 1 / num_sections) / (
                (num_sections - 2) / num_sections)
        control_name = control_surface.find('ParmContainer/Name').text
        start_length_c = 1 - float(control_surface.find('ParmContainer/SS_Control/Length_C_Start').get('Value'))
        end_length_c = 1 - float(control_surface.find('ParmContainer/SS_Control/Length_C_End').get('Value'))
        control_surfaces.append(ControlSurfaceSec(control_name, u_start, start_length_c))
        control_surfaces.append(ControlSurfaceSec(control_name, u_end, end_length_c))

    control_surfaces.sort(key=lambda ctrl: ctrl.u_loc)  # sort according to u location

    count = 0
    incidence = 0
    dihedral = 0
    cumulative_span = 0
    for XSec in geom.findall('WingGeom/XSecSurf/XSec'):
        chord = float(XSec.find('ParmContainer/XSec/Tip_Chord').get('Value'))
        root_chord = float(XSec.find('ParmContainer/XSec/Root_Chord').get('Value'))
        sweep = float(XSec.find('ParmContainer/XSec/Sweep').get('Value'))
        sweep_loc = float(XSec.find('ParmContainer/XSec/Sweep_Location').get('Value'))

        incidence_prev = incidence

        if relative_twist:
            incidence = incidence + float(XSec.find('ParmContainer/XSec/Twist').get('Value'))
        else:
            incidence = float(XSec.find('ParmContainer/XSec/Twist').get('Value'))

        if relative_dihedral:
            dihedral = dihedral + float(XSec.find('ParmContainer/XSec/Dihedral').get('Value'))
        else:
            dihedral = float(XSec.find('ParmContainer/XSec/Dihedral').get('Value'))

        if count == 0:  # openvsp's first section is always a dummy section that does not appear in the geometry
            span = 0
        else:
            span = float(XSec.find('ParmContainer/XSec/Span').get('Value'))
            cumulative_span = cumulative_span + span

        xloc_prev = xloc
        yloc_prev = yloc
        zloc_prev = zloc

        # calculate new LE coordinates
        xloc = xloc + root_chord * sweep_loc + span * math.tan(math.radians(sweep)) - chord * sweep_loc
        yloc = yloc + span * math.cos(math.radians(dihedral + xrot))
        zloc = zloc + span * math.sin(math.radians(dihedral + xrot))

        naca_airfoil = None
        afile_airfoil = None

        # set airfoil type and data
        if XSec.find('XSec/XSecCurve/XSecCurve/Type').text == '7':
            thick_chord = round(
                float(XSec.find('XSec/XSecCurve/ParmContainer/XSecCurve/ThickChord').get('Value')) * 100)
            camber = round(float(XSec.find('XSec/XSecCurve/ParmContainer/XSecCurve/Camber').get('Value')) * 100)
            if camber > 0:
                camber_loc = round(
                    float(XSec.find('XSec/XSecCurve/ParmContainer/XSecCurve/CamberLoc').get('Value')) * 10)
            else:
                camber_loc = 0
            naca_airfoil = str(camber) + str(camber_loc) + str(thick_chord)
        elif XSec.find('XSec/XSecCurve/XSecCurve/Type').text == '12':
            afile_airfoil = XSec.find('XSec/XSecCurve/FileAirfoil/AirfoilName').text + '.dat'
        else:
            naca_airfoil = '0010'

        skip_section = False
        prev_u_loc = 0

        # iterates through all control surfaces and checks if one needs to be inserted before the next section
        for control_surface in control_surfaces:
            #if control_surface.u_loc * totalSpan / 2 < cumulative_span + control_surface_tolerance and not control_surface.is_written:
            if control_surface.u_loc < count / (num_sections-2) + control_surface_tolerance and not control_surface.is_written:
                if abs(control_surface.u_loc - prev_u_loc <= control_surface_tolerance) and not count == 0:
                    avl.addControl(name, 'last', control_surface.name, control_surface.name.lower(), round(control_surface.xhinge, 3), [0, 0, 0])
                    control_surface.is_written = True

                elif abs(control_surface.u_loc - count / (num_sections-2)) <= control_surface_tolerance or count == 0:
                    # adds control surface to existing section
                    avl.addSection(name, [round(xloc, digit_precision), round(yloc, digit_precision),
                                          round(zloc, digit_precision)], round(chord, digit_precision),
                                   round(incidence, digit_precision), NACA=naca_airfoil, AFILE=afile_airfoil)
                    avl.addControl(name, 'last', control_surface.name, control_surface.name.lower(),
                                   round(control_surface.xhinge, 3), [0, 0, 0])
                    prev_u_loc = control_surface.u_loc
                    control_surface.is_written = True
                    skip_section = True
                else:
                    # defines new section for the control surface by interpolating between sections
                    #sec_u = (control_surface.u_loc * totalSpan / 2 - (cumulative_span - span)) / span
                    sec_u = (control_surface.u_loc - (count-1) / (num_sections-2)) / (1/(num_sections-2))
                    ctrl_chord = sec_u * (chord - root_chord) + root_chord
                    ctrl_ainc = sec_u * (incidence - incidence_prev) + incidence_prev
                    xloc_ctrl = sec_u * (xloc - xloc_prev) + xloc_prev
                    yloc_ctrl = sec_u * (yloc - yloc_prev) + yloc_prev
                    zloc_ctrl = sec_u * (zloc - zloc_prev) + zloc_prev

                    avl.addSection(name, [round(xloc_ctrl, digit_precision), round(yloc_ctrl, digit_precision),
                                          round(zloc_ctrl, digit_precision)], round(ctrl_chord, digit_precision),
                                   round(ctrl_ainc, digit_precision), NACA=naca_airfoil, AFILE=afile_airfoil)
                    avl.addControl(name, 'last', control_surface.name, control_surface.name.lower(),
                                   round(control_surface.xhinge, 3), [0, 0, 0])
                    prev_u_loc = control_surface.u_loc
                    control_surface.is_written = True

        if not skip_section:
            avl.addSection(name,
                           [round(xloc, digit_precision), round(yloc, digit_precision), round(zloc, digit_precision)],
                           round(chord, digit_precision),
                           round(incidence, digit_precision), NACA=naca_airfoil, AFILE=afile_airfoil)

        count = count + 1


for geom in root.findall('Vehicle/Geom'):
    if geom.find('ParmContainer/WingGeom') is not None and '*' not in geom.find('ParmContainer/Name').text: # checks for excluded wing geometry
        if geom.find('ParmContainer/Name').text == 'Wing_ref':
            avl.setRefSurf('Wing_ref')  # this is here only so that avlGeo does not throw an error
            Sref = float(geom.find('ParmContainer/WingGeom/TotalArea').get('Value'))
            Bref = float(geom.find('ParmContainer/WingGeom/TotalProjectedSpan').get('Value'))
            Cref = float(geom.find('ParmContainer/WingGeom/TotalChord').get('Value'))

            # manually set reference dimensions because I trust openvsp's method more than avlGeo
            avl.setRefs(round(Sref, digit_precision), round(Bref, digit_precision), round(Cref, digit_precision))

            totalSpan = float(geom.find('ParmContainer/WingGeom/TotalSpan').get('Value'))   # pass non-projected span
            define_surface(geom, totalSpan, Cref)
        else:
            totalSpan = float(geom.find('ParmContainer/WingGeom/TotalSpan').get('Value'))
            chord = float(geom.find('ParmContainer/WingGeom/TotalChord').get('Value'))
            define_surface(geom, totalSpan, chord)

avl.write(avl_out_path)
