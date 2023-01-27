

class avlGeo:
    def __init__(self, name, manualrefs):
        # References
        self.name = name
        self.manualrefs = manualrefs

        self.refSurf = None
        self.sref = None
        self.cref = None
        self.bref = None
        self.xref = None
        self.yref = None
        self.zref = None
        self.IYsym = None
        self.IZsym = None
        self.Zsym = None
        self.mach = None

        self.surfaces = {}

    def addSurface(self, name, Nchord, Angle, Cspace, Nspan=None, Sspace=None, Ydup=None):
        # Safety Check
        if name in list(self.surfaces.keys()):
            raise Exception('Surface with name {} already exists!'.format(name))

        self.surfaces[name] = {'Nchord': Nchord, 'Cspace': Cspace, 'YDUPLICATE': Ydup, 'Angle': Angle, 'Nspan': Nspan,
                               'Sspace': Sspace, 'Sections': []}

    def addSection(self, surf_name, xyz_le, Chord, Ainc, Nspan=None, Sspace=None, NACA=None, AFILE=None):
        # Safety Check
        self.exists(surf_name)

        # Check that either AFILE or NACA is defined.
        if NACA is None and AFILE is None:
            raise Exception('Please either provide NACA number or airfoil .dat file.')

        if NACA is not None and AFILE is not None:
            raise Exception('Please only provide either a NACA number of airfoil .dat file')

        # Check if spanwise attributes are not defined in either the surface or section.
        if self.surfaces[surf_name]['Nspan'] is None and Nspan is None:
            raise Exception('Nspan not provided in the surface definition. Please provide for this section definition')

        if self.surfaces[surf_name]['Sspace'] is None and Sspace is None:
            raise Exception('Sspace not provided in the surface definition. Please provide for this section definition')

        # Check if spanwise attributes are double defined in both the surface or section.
        if self.surfaces[surf_name]['Nspan'] is not None and Nspan is not None:
            raise Exception('Nspan defined in both the surface and section definition. Please pick one')

        if self.surfaces[surf_name]['Sspace'] is not None and Sspace is not None:
            raise Exception('Sspace defined in both the surface and section definition. Please pick one')

        if NACA is not None:
            new_section = {'Xle': xyz_le[0], 'Yle': xyz_le[1], 'Zle': xyz_le[2], 'Chord': Chord, 'Ainc': Ainc,
                           'Nspan': Nspan, 'Sspace': Sspace, 'NACA': NACA, 'Controls': {}}

        if AFILE is not None:
            new_section = {'Xle': xyz_le[0], 'Yle': xyz_le[1], 'Zle': xyz_le[2], 'Chord': Chord, 'Ainc': Ainc,
                           'Nspan': Nspan, 'Sspace': Sspace, 'AFILE': AFILE, 'Controls': {}}

        self.surfaces[surf_name]['Sections'].append(new_section)

    def addControl(self, surf_name, sec_idx, ctrl_name, control_type, x_hinge, hinge):
        # Sec_idx can also be 'last', automatically add control surface to the last section of that surface
        if isinstance(sec_idx, (str)) and sec_idx == 'last':
            sec_idx = self.numSec(surf_name) - 1
        # Safety Check
        self.exists(surf_name, sec_idx=sec_idx)

        new_control = {'type': control_type, 'x_hinge': x_hinge, 'hinge': hinge}

        self.surfaces[surf_name]['Sections'][sec_idx]['Controls'][ctrl_name] = new_control

    def spacingStr2Int(self, dict_in):
        # This will search for Cspace and Sspace in the given dictionary and translate 'cosine' into 1.0
        '''
                parameter                              spacing
        ---------                              -------

            3.0        equal         |   |   |   |   |   |   |   |   |

            2.0        sine          || |  |   |    |    |     |     |

            1.0        cosine        ||  |    |      |      |    |  ||

            0.0        equal         |   |   |   |   |   |   |   |   |

            -1.0        cosine        ||  |    |      |      |    |  ||

            -2.0       -sine          |     |     |    |    |   |  | ||

            -3.0        equal         |   |   |   |   |   |   |   |   |

            Sspace (spanwise)  :    first section        ==>       last section
            Cspace (chordwise) :    leading edge         ==>       trailing edge
            Bspace (lengthwise):    frontmost point      ==>       rearmost point
            '''

        for k in list(dict_in.keys()):
            if k == 'Cspace' or k == 'Sspace':
                if dict_in[k] == 'cosine':
                    dict_in[k] = 1
                elif dict_in[k] == 'equal':
                    dict_in[k] = 0
                elif dict_in[k] == 'sine':
                    dict_in[k] = 2
                elif dict_in[k] == '-sine':
                    dict_in[k] = -2
        return dict_in

    def write(self, fname):
        if self.refSurf is None:
            raise Exception(
                'Reference surface not set. Please run setRefSurf() and provide the name of the surface to use for reference values.')
        self.updateRefs()

        with open(fname, 'w') as fh:
            fh.write(self.name)
            fh.write('\n')
            fh.write('#MACH\n')
            fh.write(str(self.mach))
            fh.write('\n')
            # Write symetery
            fh.write('#IYsym 	IZsym 	Zsym 	Vehicle Symmetry\n')
            fh.write(str(self.IYsym) + ' ')
            fh.write(str(self.IZsym) + ' ')
            fh.write(str(self.Zsym) + ' ')
            fh.write('\n')
            # write reference values
            fh.write('#Sref 	Cref 	Bref 	Reference Area and Lengths\n')
            fh.write(str(self.sref) + ' ')
            fh.write(str(self.cref) + ' ')
            fh.write(str(self.bref) + ' ')
            fh.write('\n')
            # write cg location
            fh.write('#Xref 	Yref 	Zref 	Center of Gravity Location\n')
            fh.write(str(self.xref) + ' ')
            fh.write(str(self.yref) + ' ')
            fh.write(str(self.zref) + ' ')
            fh.write('\n')
            for surf_name in list(self.surfaces.keys()):
                fh.write('\n\n#===============================\n')
                fh.write('SURFACE\n')
                fh.write(surf_name + '\n')

                # Write SURFACE attribute names
                fh.write('#Nchord Cspace Nspan Sspace\n')
                # Write SURFACE attribute values
                self.surfaces[surf_name] = self.spacingStr2Int(self.surfaces[surf_name])
                fh.write(str(self.surfaces[surf_name]['Nchord']) + ' ')
                fh.write(str(self.surfaces[surf_name]['Cspace']) + ' ')
                if self.surfaces[surf_name]['Nspan'] is not None:
                    fh.write(str(self.surfaces[surf_name]['Nspan']) + ' ')
                if self.surfaces[surf_name]['Sspace'] is not None:
                    fh.write(str(self.surfaces[surf_name]['Sspace']) + ' ')
                fh.write('\n')
                if self.surfaces[surf_name]['YDUPLICATE'] is not None:
                    fh.write('YDUPLICATE\n')
                    fh.write(str(self.surfaces[surf_name]['YDUPLICATE']))
                fh.write('\n')
                fh.write('ANGLE\n')
                fh.write(str(self.surfaces[surf_name]['Angle']))
                fh.write('\n')

                for section in self.surfaces[surf_name]['Sections']:
                    sec = self.spacingStr2Int(section)
                    fh.write('#------------------------\n')
                    fh.write('SECTION\n')
                    # Write SECTION attribute names
                    fh.write('#Xle Yle Zle Chord Ainc Nspan Sspace\n')
                    # Write SECTION attributes
                    fh.write(str(sec['Xle']) + ' ')
                    fh.write(str(sec['Yle']) + ' ')
                    fh.write(str(sec['Zle']) + ' ')
                    fh.write(str(sec['Chord']) + ' ')
                    fh.write(str(sec['Ainc']) + ' ')
                    if sec['Nspan'] is None:
                        sec['Nspan'] = 0
                    fh.write(str(sec['Nspan']) + ' ')
                    if sec['Sspace'] is None:
                        sec['Sspace'] = 0
                    fh.write(str(sec['Sspace']) + ' ')
                    fh.write('\n')
                    if 'AFILE' in list(sec.keys()):
                        fh.write('AFILE\n')
                        fh.write(sec['AFILE'])
                    elif 'NACA' in list(sec.keys()):
                        fh.write('NACA\n')
                        fh.write(sec['NACA'])
                    fh.write('\n')
                    # Write control surfaces
                    for ctrl in list(sec['Controls'].keys()):
                        fh.write('CONTROL\n')
                        fh.write('#NAME GAIN Xhinge XHVEC YHVEC ZHVEC SGNDUP\n')
                        fh.write(ctrl + ' ')
                        if sec['Controls'][ctrl]['type'] == 'elevator' or sec['Controls'][ctrl]['type'] == 'flap':
                            fh.write('1 ')
                        elif sec['Controls'][ctrl]['type'] == 'aileron' or sec['Controls'][ctrl]['type'] == 'rudder':
                            fh.write('-1 ')
                        fh.write(str(sec['Controls'][ctrl]['x_hinge']) + ' ')
                        fh.write(str(sec['Controls'][ctrl]['hinge'][0]) + ' ')
                        fh.write(str(sec['Controls'][ctrl]['hinge'][1]) + ' ')
                        fh.write(str(sec['Controls'][ctrl]['hinge'][2]) + ' ')
                        if sec['Controls'][ctrl]['type'] == 'elevator' or sec['Controls'][ctrl]['type'] == 'flap':
                            fh.write('1 ')
                        elif sec['Controls'][ctrl]['type'] == 'aileron' or sec['Controls'][ctrl]['type'] == 'rudder':
                            fh.write('-1')
                        fh.write('\n')

    def numSec(self, surf_name):
        return len(self.surfaces[surf_name]['Sections'])

    def modSurface(self, surf_name, attr, new_value):
        self.exists(surf_name, attr=attr)
        self.surfaces[surf_name][attr] = new_value

    def modSection(self, surf_name, sec_idx, attr, new_value):
        self.exists(surf_name, sec_idx=sec_idx, attr=attr)
        self.surfaces[surf_name]['Sections'][sec_idx][attr] = new_value

    def modCtrl(self, surf_name, sec_idx, ctrl_name, attr, new_value):
        self.exists(surf_name, sec_idx=sec_idx, ctrl_name=ctrl_name, attr=attr)
        self.surfaces[surf_name]['Sections'][sec_idx]['Controls'][ctrl_name][attr] = new_value

    def exists(self, surf_name, sec_idx=None, ctrl_name=None, attr=None):
        if surf_name not in list(self.surfaces.keys()):
            raise Exception('Surface {} does not exist!'.format(surf_name))
        if sec_idx is None:  # Not looking for section attribute
            if attr is not None and attr not in list(self.surfaces[surf_name].keys()):
                raise Exception('Attribute {} does not exist for surface {}'.format(attr, surf_name))
        else:  # Looking for section index
            if sec_idx >= self.numSec(surf_name):
                raise Exception(
                    'Surface has {} sections, but attempted to reference {}th section! Note: Sections index from 0.'.format(
                        self.numSec(surf_name), sec_idx + 1))

            if ctrl_name is None:  # Not looking for a control attribute
                if attr is not None and attr not in list(self.surfaces[surf_name]['Sections'][sec_idx].keys()):
                    raise Exception(
                        'Attribute {} does not exist for section {} of surface {} '.format(attr, sec_idx, surf_name))

            else:
                if ctrl_name not in list(self.surfaces[surf_name]['Sections'][sec_idx]['Controls'].keys()):
                    raise Exception(
                        'Control named {} does not exist in section {} of surface {}'.format(ctrl_name, sec_idx,
                                                                                             surf_name))

                if attr is not None and attr not in list(
                        self.surfaces[surf_name]['Sections'][sec_idx]['Controls'][ctrl_name].keys()):
                    raise Exception(
                        'Attribute {} does not exist for control name {} of section {} of surface {}'.format(attr,
                                                                                                             ctrl_name,
                                                                                                             sec_idx,
                                                                                                             surf_name))

    def setRefSurf(self, surf_name):
        self.refSurf = surf_name

    def updateRefs(self):
        if not self.manualrefs:
            secs = self.surfaces[self.refSurf]['Sections']
            dys = [abs(secs[i + 1]['Yle'] - secs[i]['Yle']) for i in range(self.numSec(self.refSurf) - 1)]
            cs = [s['Chord'] for s in secs]
            areas = [(cs[i] + cs[i + 1]) / 2 * dys[i] for i in range(self.numSec(self.refSurf) - 1)]
            self.bref = sum(dys)
            self.sref = sum(areas)
            self.cref = sum(areas) / sum(dys)
            if self.surfaces[self.refSurf]['YDUPLICATE'] is not None:
                self.bref = self.bref * 2
                self.sref = self.sref * 2

    def setRefs(self, sref, bref, cref):
        self.sref = sref
        self.bref = bref
        self.cref = cref




