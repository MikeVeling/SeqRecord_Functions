###############################################################################
#                            Import dependencies                              #
###############################################################################
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition, SeqFeature, FeatureLocation
from collections import OrderedDict
from warnings import warn
import sys, warnings, copy, re
###############################################################################
#                            Internal Functions                               #
###############################################################################
def feature_warning():
    '''
    Warns the user that they are altering sequence in a way that overlaps with
    existing features. This would normally delete the feature as it is no 
    longer the same, but I have optionally allowed users to maintain such
    features
    '''
    warn('you are trying to edit sequence inside a feature, I have '+
         'maintained this feature annotation, but this annotation may no '+
         'longer be accurate as sequence was removed from it')
def format_colors_for_benchling(SeqRecord_obj):
    '''
    This function is largely internal, it only reformats the color 
    specification so that it can be understood by both Benchling and SnapGene.
    '''
    SeqRecord_obj_ouptut=copy.deepcopy(SeqRecord_obj)
    for feature in SeqRecord_obj_ouptut.features:
        color_in_note=False
        color='N/A'
        if 'note' in feature.qualifiers:
            for note in feature.qualifiers['note']:
                if 'color: ' in note:
                    color_in_note=True
                    if color == 'N/A':
                        color=note.split('color: ')[1]
                    elif note.split('color: ')[1] != color:
                        warn('one of your features has multiple colors')
        if 'ApEinfo_fwdcolor' in feature.qualifiers:
            if color == 'N/A':
                color=feature.qualifiers['ApEinfo_fwdcolor'][0]
            elif feature.qualifiers['ApEinfo_fwdcolor'][0] != color:
                warn('one of your features has multiple colors')
        if 'ApEinfo_revcolor' in feature.qualifiers:
            if color == 'N/A':
                color=feature.qualifiers['ApEinfo_revcolor'][0]
            elif feature.qualifiers['ApEinfo_revcolor'][0] != color:
                warn('one of your features has multiple colors')
        if color != 'N/A':
            if not color_in_note:
                if 'note' not in feature.qualifiers:
                    feature.qualifiers['note']=[]
                feature.qualifiers['note'].append('color: '+color)
            feature.qualifiers['ApEinfo_fwdcolor']=[color]
            feature.qualifiers['ApEinfo_revcolor']=[color]
    return SeqRecord_obj_ouptut
###############################################################################
#                               IO Functions                                  #
###############################################################################
def open_SeqRecord_obj(path):
    '''
    Open a .gbk file and return a SeqRecord object
    '''
    return format_colors_for_benchling(next(SeqIO.parse(path, format='gb')))

def save_SeqRecord_obj(SeqRecord_obj,
                 path):
    '''   
    Save a SeqRecord object to a .gbk file
    '''
    SeqIO.write(SeqRecord_obj,path,format='gb')
###############################################################################
#                              Base Functions                                 #
###############################################################################
def find_sequence(SeqRecord_obj,
                  search_seq,
                  forward_only=False,
                  reverse_only=False):
    '''
    Returns all locations of your exact search sequence. If directionality is
    important, use the toggles.
    '''
    if search_seq == '':
        return []
    forward_locs=[m.start() for m in re.finditer('(?='+search_seq.upper()+')',
                                               str(SeqRecord_obj.seq).upper())]
    rc_search_seq=str(Seq(search_seq,
                      SeqRecord_obj.seq.alphabet).reverse_complement()).upper()
    reverse_locs=[m.start() for m in re.finditer('(?='+rc_search_seq+')',
                                               str(SeqRecord_obj.seq).upper())]
    all_locs=[]
    for loc in forward_locs:
        if not reverse_only:
            all_locs.append([loc,loc+len(search_seq),False])
    for loc in reverse_locs:
        if not forward_only:
            all_locs.append([loc,loc+len(search_seq),True])
    return all_locs

def remove_sequence(SeqRecord_obj,
                    start_input,
                    end_input,
                    maintain_ambiguious_features=True,
                    supress_warnings=False):
    '''
    Will remove a sequence based on a particular region within the sequence.
    This function takes a SeqRecord_obj in the form of a SeqRecord and will 
    remove the sequence between the start and end input location. By default, 
    it will not remove features that are edited, but can remove such features
    if maintain_ambiguious_features=False. It will warn users about editing 
    such features unless supress_warnings=True. This is a base function that is
    abstracted upon in other functions.
    
    Note: This function is strand agnostic. Therefore start must be before 
    end
    '''
    if end_input < start_input:
        start=end_input
        end=start_input
    else:
        start=start_input
        end=end_input
    SeqRecord_obj_intermediate=copy.deepcopy(SeqRecord_obj)
    SeqRecord_obj_ouptut=(SeqRecord_obj_intermediate[0:start]+
                          SeqRecord_obj_intermediate[end:])
    SeqRecord_obj_ouptut.dbxrefs = SeqRecord_obj_intermediate.dbxrefs
    SeqRecord_obj_ouptut.annotations = SeqRecord_obj_intermediate.annotations
    if maintain_ambiguious_features:
        for feature_input in SeqRecord_obj_intermediate.features:
            feature=copy.deepcopy(feature_input)
            feature_start=int(feature.location.start)
            feature_end=int(feature.location.end)
            if (feature_start <= start and
               feature_end >= end) and not (
               feature_start == start and
               feature_end == end):
                if not supress_warnings:
                    feature_warning()
                feature.location._end=ExactPosition(feature_end-(end-start))
                SeqRecord_obj_ouptut.features.append(feature)
            elif feature_start == start and (
                  feature_end == end):
                pass
            elif (feature_start < start and
                 feature_end > start and
                 feature_end <= end):                                          # End spanning
                feature.location._end=ExactPosition(start)
                if not supress_warnings:
                    feature_warning()
                SeqRecord_obj_ouptut.features.append(feature)
            elif (start <= feature_start and
                  end >= feature_end):                                         # Fully spanning case
                pass
            elif (feature_start >= start and
                 feature_start < end and
                 feature_end > end):                                           # start spanning
                if not supress_warnings:
                    feature_warning()
                feature.location._end=ExactPosition(feature_end-(end-start))
                feature.location._start=ExactPosition(start)
                SeqRecord_obj_ouptut.features.append(feature)
            elif (feature_start >= end and
                 feature_end >= end):                                          # After case
                 pass
            elif (feature_start < start and
                 feature_end <= start):                                        # Fully Before
                pass
            else:
                print('wut happened? you are trying to edit a sequence in a'+
                      ' way that overlaps with features in a way I did not'+
                      ' expect, check the code')
                print(feature_input)
                sys.exit(1)
    return SeqRecord_obj_ouptut

def add_sequence(SeqRecord_obj,
                 start,
                 sequence,
                 make_feature=True,
                 feature_name='N/A',
                 feature_type='CDS',
                 feature_qualifiers=OrderedDict(),
                 feature_color='#FF0000',
                 strand=1,
                 maintain_ambiguious_features=True,
                 supress_warnings=False):
    '''
    Add a sequence to a specified start location within a SeqRecord_obj as a
    SeqRecord object. It can also make that newly added sequence a feature (as
    default). If a feature is made, passthrough arguments can be used to
    specify feature name, feature type, and feature color. the 
    feature_qualifiers are a set of other data that could add extra information
    that different visualizers use to interpret the feature, but is unnecessary
    in most cases. The same consideration about editing sequence within other
    features is here (see above feature_warning function). This is a base
    function that is abstracted upon with other functions.
    
    Note: this function is not strand agnostic. by default it adds your 
    sequence in the forward direction, but if strand=-1 the sequence will be 
    added in reverse at the same start location.
    '''
    SeqRecord_obj_intermediate=copy.deepcopy(SeqRecord_obj)
    if isinstance(sequence, str):
        sequence_useful=Seq(sequence,
                            alphabet=SeqRecord_obj_intermediate.seq.alphabet)
    else:
        sequence_useful=sequence
    if strand == 1:
        SeqRecord_obj_ouptut=(SeqRecord_obj_intermediate[0:start]+
                        sequence_useful+
                        SeqRecord_obj_intermediate[start:])
    elif strand == -1:
        SeqRecord_obj_ouptut=(SeqRecord_obj_intermediate[0:start]+
                        sequence_useful.reverse_complement()+
                        SeqRecord_obj_intermediate[start:])
    SeqRecord_obj_ouptut.dbxrefs = SeqRecord_obj_intermediate.dbxrefs
    SeqRecord_obj_ouptut.annotations = SeqRecord_obj_intermediate.annotations
    if maintain_ambiguious_features:
        for feature_input in SeqRecord_obj_intermediate.features:
            feature=copy.deepcopy(feature_input)
            feature_start=int(feature.location.start)
            feature_end=int(feature.location.end)
            if (feature_start < start and
                feature_end > start):
                if not supress_warnings:
                    feature_warning()
                feature.location._end=ExactPosition(feature_end+
                                                    len(sequence_useful))
                SeqRecord_obj_ouptut.features.append(feature)
    if make_feature:
        SeqRecord_obj_ouptut=add_feature(SeqRecord_obj_ouptut,
                                   start,
                                   start+len(sequence_useful),
                                   feature_name=feature_name,
                                   feature_type=feature_type,
                                   feature_qualifiers_input=feature_qualifiers,
                                   feature_color=feature_color,
                                   strand=strand,
                                   supress_warnings=supress_warnings)
    return SeqRecord_obj_ouptut
###############################################################################
#                            Feature Functions                                #
###############################################################################
def identify_feature(SeqRecord_obj,
                     feature):
    '''
    finds a list of things that could be your feature. If you give it a string,
    it will find all features with that name. If you give it an individual
    feature, it will make sure that feature is in your SeqRecord_obj and return
    a list containing just that feature. If it is a list of feature objects 
    already, it will just return that same list.
    '''
    if isinstance(feature,str):
        target_features=[]
        for feature_obj in SeqRecord_obj.features:
            if 'label' in feature_obj.qualifiers:
                feature_name=feature_obj.qualifiers['label']
                if feature in feature_name:
                    target_features.append(feature_obj)
    elif isinstance(feature, SeqFeature):
        if feature not in SeqRecord_obj.features:
            raise Exception ('Your target feature is not in your SeqRecord '+
                             'object')
        target_features=[feature]
    elif isinstance(feature, list):
        for feature_test in feature:
            if feature_test not in SeqRecord_obj.features:
                raise Exception ('Your target feature is not in your '+
                                 'SeqRecord object')
        target_features=feature
    else:
        raise Exception(feature+' can not be parsed into feature objects')
    return target_features

def add_translation_to_feature(SeqRecord_obj,
                               Feature_obj,
                               supress_warnings=False):
    '''
    Will add a translation to a feature object in a format understood by 
    snapgene. Note this is not done in place, so make sure to add the feature
    back to your SeqRecord_obj
    '''
    Feature_obj_output=copy.deepcopy(Feature_obj)
    feature_qualifiers=Feature_obj_output.qualifiers
    if 'translation' not in feature_qualifiers:
        start=int(Feature_obj_output.location.start)
        end=int(Feature_obj_output.location.end)
        assert end >= start
        strand=Feature_obj_output.location.strand
        if supress_warnings:
            with warnings.catch_warnings():
                if strand == 1:
                    feature_qualifiers['translation']=str(
                                      SeqRecord_obj[start:end].translate().seq)
                elif strand == -1:
                    feature_qualifiers['translation']=str(
                 SeqRecord_obj[start:end].reverse_complement().translate().seq)
                else:
                    raise Exception('Strand not understood')
        else:
            if strand == 1:
                feature_qualifiers['translation']=str(
                                      SeqRecord_obj[start:end].translate().seq)
            elif strand == -1:
                feature_qualifiers['translation']=str(
                 SeqRecord_obj[start:end].reverse_complement().translate().seq)
            else:
                raise Exception('Strand not understood')
    return Feature_obj_output

def add_feature(SeqRecord_obj,
                start,
                end,
                feature_name='N/A',
                feature_type='CDS',
                feature_qualifiers_input=OrderedDict(),
                feature_color='#FF0000',
                strand=1,
                supress_warnings=False):
    '''
    Adds a feature annotation to a SeqRecord_obj which is a SeqRecord object. 
    This feature is added between the start and end location specified. Feature
    name, type, and color can all be specified. the feature_qualifiers can be
    used to pass an ordered dictionary containing extra information about the 
    feature directly to the feature. This extra information could be used by 
    other programs to glen more information about that feature, but is largely
    unnecessary. If the feature_type is a coding sequence ('CDS'), a 
    translation is included. This is a base function that is abstracted upon
    in other functions here.
    
    Note: this function does not alter sequence, it simply adds a feature 
    annotation to existing sequence.
    '''
    feature_qualifiers=copy.deepcopy(feature_qualifiers_input)
    assert start < end
    if 'label' not in feature_qualifiers:
        feature_qualifiers['label']=[]
    if 'note' not in feature_qualifiers:
        feature_qualifiers['note']=[]
    if 'ApEinfo_fwdcolor' not in feature_qualifiers:
        feature_qualifiers['ApEinfo_fwdcolor']=[]
    if 'ApEinfo_revcolor' not in feature_qualifiers:
        feature_qualifiers['ApEinfo_revcolor']=[]
    feature_qualifiers['label'].append(feature_name)
    feature_qualifiers['note'].append('color: '+feature_color)
    feature_qualifiers['ApEinfo_fwdcolor'].append(feature_color)
    feature_qualifiers['ApEinfo_revcolor'].append(feature_color)
    new_feature=SeqFeature(location=FeatureLocation(start,end),
                           type=feature_type,
                           strand=strand,
                           qualifiers=feature_qualifiers)
    if feature_type == 'CDS':
        new_feature=add_translation_to_feature(SeqRecord_obj,
                                               new_feature,
                                               supress_warnings=(
                                                             supress_warnings))
    SeqRecord_obj_output=copy.deepcopy(SeqRecord_obj)
    SeqRecord_obj_output.features.append(new_feature)
    return SeqRecord_obj_output

def remove_feature(SeqRecord_obj,
                   feature):
    '''
    Removes a feature by name or by direct passing of the feature object.
    '''
    SeqRecord_obj_output=copy.deepcopy(SeqRecord_obj)
    target_features=identify_feature(SeqRecord_obj_output,feature)
    if len(target_features) == 0:
        warn('no features have been removed as we could not identify your '+
             'feature')
    elif len(target_features) > 1:
        warn('multiple features have the same name, I will delete all of them')
    for feature_obj in target_features:
        SeqRecord_obj_output.features.remove(feature_obj)
    return SeqRecord_obj_output
###############################################################################
#                          Abstracted Functions                               #
###############################################################################
def add_feature_from_sequence(SeqRecord_obj,
                              search_seq,
                              feature_name='N/A',
                              feature_type='CDS',
                              feature_qualifiers_input=OrderedDict(),
                              feature_color='#FF0000',
                              forward_only=False,
                              reverse_only=False,
                              supress_warnings=False):
    '''
    adds feature using a search sequence. If multiple copies of the feature are
    present, it will add multiple copies. This is an abstraction of the
    add_feature function. 
    '''
    SeqRecord_obj_output=copy.deepcopy(SeqRecord_obj)
    locations=find_sequence(SeqRecord_obj_output,
                            search_seq,
                            forward_only=forward_only,
                            reverse_only=reverse_only)
    if not supress_warnings:
        if len(locations) == 0:
            warn('could not find the sequence of the feature you are trying '+
                 'to add. No features were added')
        elif len(locations) > 1:
            warn('Found many copies of your sequence, adding features to all '+
                 'of them')
    for location in locations:
        start=location[0]
        end=location[1]
        if location[2]:
            strand=-1
        else:
            strand=1
        SeqRecord_obj_output=add_feature(SeqRecord_obj_output,
                                   start,
                                   end,
                                   feature_name=feature_name,
                                   feature_type=feature_type,
                             feature_qualifiers_input=feature_qualifiers_input,
                                   feature_color=feature_color,
                                   strand=strand,
                                   supress_warnings=supress_warnings)
    return SeqRecord_obj_output
def delete_seq_by_features(SeqRecord_obj,
                           feature_input,
                           maintain_ambiguious_features=True,
                           supress_warnings=False):
    '''
    Remove a feature and the underlying sequence. This is an abstraction of
    the remove_sequence.
    '''
    SeqRecord_obj_output=copy.deepcopy(SeqRecord_obj)
    target_features=identify_feature(SeqRecord_obj_output,feature_input)
    if not supress_warnings:
        if len(target_features) == 0:
            warn('could not find the feature you are trying '+
                 'to delete. No features were deleted')
        elif len(target_features) > 1:
            warn('Found many copies of your feature, deleting all of them')
    while len(target_features) > 0:
        feature=target_features[0]
        start=int(feature.location.start)
        end=int(feature.location.end)
        SeqRecord_obj_output=remove_sequence(SeqRecord_obj_output,
                                       start,
                                       end,
                                       maintain_ambiguious_features=(
                                                 maintain_ambiguious_features),
                                       supress_warnings=supress_warnings)
        target_features=identify_feature(SeqRecord_obj_output,feature_input)
    return SeqRecord_obj_output
def replace_sequence_by_locus(SeqRecord_obj,
                              start,
                              end,
                              sequence,
                              strand=1,
                              make_feature=False,
                              feature_name='N/A',
                              feature_type='CDS',
                              feature_qualifiers=OrderedDict(),
                              feature_color='#FF0000',
                              maintain_ambiguious_features=True,
                              supress_warnings=False):
    '''
    Will replace a sequence specified by the location with a new input sequence
    This is an abstraction of the remove_sequence and add_sequence functions.
    '''
    SeqRecord_obj_intermediate=copy.deepcopy(SeqRecord_obj)
    SeqRecord_obj_intermediate=remove_sequence(SeqRecord_obj_intermediate,
                                         start,
                                         end,
                                         maintain_ambiguious_features=(
                                                 maintain_ambiguious_features),
                                         supress_warnings=supress_warnings)
    SeqRecord_obj_output=add_sequence(SeqRecord_obj_intermediate,
                                start,
                                sequence,
                                make_feature=make_feature,
                                feature_name=feature_name,
                                feature_type=feature_type,
                                feature_qualifiers=feature_qualifiers,
                                feature_color=feature_color,
                                strand=strand,
                                maintain_ambiguious_features=(
                                                 maintain_ambiguious_features),
                                supress_warnings=supress_warnings)
    return SeqRecord_obj_output
def replace_sequence_by_feature(SeqRecord_obj,
                                target_feature,
                                sequence,
                                make_feature=False,
                                feature_name='N/A',
                                feature_type='CDS',
                                feature_qualifiers=OrderedDict(),
                                feature_color='#FF0000',
                                maintain_ambiguious_features=True,
                                supress_warnings=False):
    '''
    Will take a target feature and replace all instances of that target feature
    with a new sequence. this is an abstraction of replace_sequence_by_locus.
    '''
    SeqRecord_obj_intermediate=copy.deepcopy(SeqRecord_obj)
    target_features=identify_feature(SeqRecord_obj_intermediate,target_feature)
    if len(target_features) == 0:
        warn('no features have been replaced as we could not identify your '+
             'feature')
    elif len(target_features) > 1:
        warn('multiple features have the same name, I will replace all of '+
             'them')
    for feature in target_features:
        start=int(feature.location.start)
        end=int(feature.location.end)
        strand=feature.location.strand
        SeqRecord_obj_output=replace_sequence_by_locus(
                                                 SeqRecord_obj_intermediate,
                                                 start,
                                                 end,
                                                 sequence,
                                                 strand=strand,
                                                 make_feature=make_feature,
                                                 feature_name=feature_name,
                                                 feature_type=feature_type,
                                                 feature_qualifiers=(
                                                           feature_qualifiers),
                                                 feature_color=feature_color,
                                                 maintain_ambiguious_features=(
                                                 maintain_ambiguious_features),
                                                 supress_warnings=(
                                                             supress_warnings))
    return SeqRecord_obj_output
def translate_CDS(SeqRecord_obj,
                  supress_warnings=False):
    '''
    Will iterate through a SeqRecord_obj to add translations to all the 
    features with type "CDS"
    '''
    SeqRecord_obj_output=copy.deepcopy(SeqRecord_obj)
    for feature in SeqRecord_obj_output.features:
        if feature.type == 'CDS':
            new_feature=add_translation_to_feature(SeqRecord_obj,
                                                   feature,
                                                   supress_warnings=(
                                                             supress_warnings))
            SeqRecord_obj_output.features.remove(feature)
            SeqRecord_obj_output.features.append(new_feature)
    return SeqRecord_obj_output
###############################################################################
#                           Tests the functions                               #
###############################################################################
def test():
    '''
    This test function ensures that all various  functions work properly. it
    must be run in a folder with test_in.gbk file.
    
    Note: The expected behavior includes warning about editing features in an
    ambiguous way. Their presence with the production of new gbk files is a 
    successful run.
    '''
    SeqRecord_obj=open_SeqRecord_obj('test_in.gbk')
    SeqRecord_obj_cut=remove_sequence(SeqRecord_obj,1000,1500)
    save_SeqRecord_obj(SeqRecord_obj_cut,'test_out_cut.gbk')
    SeqRecord_obj_added=add_sequence(SeqRecord_obj,1000,'ATC')
    save_SeqRecord_obj(SeqRecord_obj_added,'test_out_added.gbk')
    SeqRecord_obj_cut_and_added=add_sequence(SeqRecord_obj_cut,1000,'ATC')
    save_SeqRecord_obj(SeqRecord_obj_cut_and_added,
                       'test_out_cut_and_added.gbk')
    SeqRecord_obj_mod_atgs=add_feature_from_sequence(
                                               SeqRecord_obj_cut_and_added,
                                               'ATG',
                                               feature_name='Start')
    save_SeqRecord_obj(SeqRecord_obj_mod_atgs,'test_moded_ATGs.gbk')
    SeqRecord_obj_remove_atgs=delete_seq_by_features(SeqRecord_obj_mod_atgs,
                                                     'Start')
    save_SeqRecord_obj(SeqRecord_obj_remove_atgs,'test_removed_ATGs.gbk')
    replaced_GFP=replace_sequence_by_feature(SeqRecord_obj_remove_atgs,
                                             'replace GFP',
                                             'CATCATCAT',
                                             make_feature=True,
                                             feature_name='kitty kat his tag')
    save_SeqRecord_obj(replaced_GFP,'replaced.gbk')
###############################################################################
#                                     Fin                                     #
###############################################################################
