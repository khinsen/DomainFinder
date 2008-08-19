import MMTK

class ParseError(Exception):
    pass

def parseRigidBodyDefinition(d, protein):
    rigid_body_groups = []
    for rigid_body_definition in d.split(','):
        rb = MMTK.Collection()
        for segment_definition in rigid_body_definition.split('+'):

            segment = segment_definition.split(':')
            if len(segment) == 2:
                chain = int(segment[0])-1
                del segment[0]
            elif len(segment) == 1:
                chain = 0
            else:
                raise ParseError("Invalid segment specification: %s\n"
                                 % segment_definition)
            if chain < 0 or chain > len(protein)+1:
                raise ParseError(("Invalid chain number: %d\n" % (chain+1)) + 
                                 ("Protein has %d chains" % len(protein)))

            residues = [int(v) for v in segment[0].split('-')]
            if len(residues) < 1 or len(residues) > 2 \
               or residues[0] > residues[1] or residues[0] < 1 \
               or residues[1] < 1:
                    raise ParseError("Invalid segment specification: %s\n"
                                     % segment_definition)
            for resnum in residues:
                if resnum-1 > len(protein[chain]):
                    raise ParseError(("Invalid residue number: %d\n"
                                      % resnum) + 
                                     ("Chain has %d residues"
                                      % len(protein[chain])))
            if len(residues) == 1:
                rb.addObject(protein[chain][residues[0]-1])
            else:
                rb.addObject(protein[chain][residues[0]-1:residues[1]])
        rigid_body_groups.append(rb)
    return rigid_body_groups

