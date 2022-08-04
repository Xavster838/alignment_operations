import pysam
from re import findall
from itertools import starmap
from cigar_opt_class import CigarOperation
from cigar_opt_class import cigar_dict

def get_cigarTuple_alignment_score(cigar_tuples):
    '''given a tuple of cigar operations (opt, l) return the score of these operations.
    '''
    cigar_opt = CigarOperation()
    scores = list( starmap(cigar_opt.operation, cigar_tuples) )
    return sum(scores)

def get_cigar_alignment_score(cigar_str):
    '''given a cigar string, process each operation and get the alignment score'''
    cigar_tuples = zip(findall('[^\d+]', cigar_str), findall('\d+', cigar_str) ) #split into operations and length and zip into tuples.
    return get_cigarTuple_alignment_score(cigar_tuples)

def get_AlignedSegment_alignment_score(aln_seg):
    '''given a pysam AlignedSegment object, compute the alignment score.'''
    assert isinstance(aln_seg, pysam.AlignmentSegment) , "Attempting to process variable that is not a pysam AlignmentSegment"
    return get_cigar_alignment_score(aln_seg.cigarstring)

def get_cigarLoc_from_rCoord( aln_seg , ref_coord ):
    '''given a pysam AlignmentSegment and a ref_coordinate, return the location on the cigarstring corresponding to that coordinate.
       @input: aln_seg : pysam AlignmentSegment.
       @input: ref_coordinate : integer value of the coordinate on the reference.
       @output: tuple of: 1. int index of cigar operation in cigarString tuple. 2. distance to move up length to get to specific coordinate
    '''
    conRef_opts = [ 'M', 'D', 'N', 'E', 'X' ]
    r_loc = aln_seg.reference_start
    assert ref_coord in range(aln_seg.reference_start, aln_seg.reference_end+1) , f"get_cigarLoc_from_rCoord Error: passed reference_coordinate {ref_coord} not in reference range: {aln_seg.reference_start}-{aln_seg.reference_end}"
    for i, cur_opt in enumerate( aln_seg.cigartuples ):
        opt , l = cur_opt[0] , cur_opt[1]
        if(cigar_dict[opt] in conRef_opts):
            l = min(l, ref_coord - r_loc)
            r_loc += l
        if(r_loc == ref_coord):
            return (i, l)
    raise Exception(f"Got to end of cigar tuple without reaching ref_coordinate {seg.reference_name} : {ref_coord}")

def get_qCoord_from_cigarLoc( aln_seg , cigar_loc):
    '''given a pysam AlignmentSegment and a ref_coordinate, return the coordinate of the query. For hard or softmasked cigarstrings.
       @input : aln_seg : see func get_cigarLoc_from_rCoord
       @input: cigar_loc : tuple of: 1. int index of cigar operation in cigarString tuple. 
                                     2. distance to move up length to get to specific coordinate
       @output : q_coordinate : coordinate in query corresponding to cigar location in the alignment.
    '''
    q_step = 0  #how far query coordinate will be from query start
    #assert aln_seg.cigarstring.find("S") == -1 , "Error: Bam needs to be hard clipped, not soft clipped." 
    clip = [0,0]
    clip[0] = aln_seg.cigar[0][1] if cigar_dict[ aln_seg.cigar[0][0] ] in ['H' , 'S'] else 0
    clip[1] = aln_seg.cigar[-1][1] if cigar_dict[ aln_seg.cigar[-1][0] ] in ['H' , 'S'] else 0
    conQuery = ['M', 'I', 'E', 'X'] #removing soft clip and treating as hard clip.
    for i in range(0,cigar_loc[0]): #run through all operations up 
        opt , l = aln_seg.cigartuples[i]
        if( cigar_dict[opt] in conQuery ):
            q_step += l
    q_step += cigar_loc[1]
    if(aln_seg.is_reverse):
        q_coord = clip[1] + aln_seg.query_alignment_length - q_step
    else:
        q_coord = clip[0] + q_step
    assert q_coord in range(clip[0] , clip[0] + aln_seg.query_alignment_length + clip[1] ) , f"problem with arithmetic. q_coord not in than query range : q_coord : {q_coord} ; alignment range : [ {clip[0]} , {clip[0] + aln_seg.query_alignment_length + clip[1]} ]"
    return q_coord


def get_qCoord_from_rCoord(aln_seg , ref_coord):
    ''' given a pysamAlignmentSegment and ref_coordinate, returen the query coordinate cooresponding to that location.
        @input: aln_seg, ref_coord (see function get_cigarLoc_from_rCoord)
        @output : q_coordinate : coordinate in query corresponding to ref_location in the alignment.
    '''
    cigar_loc = get_cigarLoc_from_rCoord( aln_seg , ref_coord )
    return get_qCoord_from_cigarLoc( aln_seg , cigar_loc )

def subset_cigar(aln_seg ,ref_start , ref_end):
    '''given reference coordinates, subset a cigarstring and return the subset cigar
       @input : aln_seg (pysam alignment segment) , ref_start (integer) , ref_end (integer)
       @output : cigar tuples list of subset
    '''
    opt_1 = get_cigarLoc_from_rCoord( aln_seg , ref_start  )
    opt_2 = get_cigarLoc_from_rCoord( aln_seg , ref_end  )
    subset_tuples = aln_seg.cigartuples[ opt_1[0] : opt_2[0] + 1 ]    #add 1 at the end because range doens't include last element.
    subset_tuples[0] , subset_tuples[-1] = (subset_tuples[0][0] , opt_1[1] )  , ( subset_tuples[-1][0] , opt_2[1] )  #change first and last operation lengths to be what get_cigarLoc output ; also need to reassign completely because can't reassign tuples
    return subset_tuples




