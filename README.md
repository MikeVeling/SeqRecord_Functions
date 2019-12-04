# SeqRecord-Functions
This repository contains my own functions built on top of the Biopython SeqRecord objects. These functions allow for easy manipulation of these objects within python. These objects can be easily imported and exported from the standard GenBank file format for annotated sequences. This format is widely used by synthetic biologists for gene design and manipulation and integrates nicely into standard visualization tools like SnapGene and Benchling.

**Note**: These functions do not consider things like cloning workflow. If you want tools to generate cloning workflows, please check out [pydna](https://github.com/BjornFJohansson/pydna)

# Installation 
These functions were written in Python 3.7.3 but should be largely compatible with any and all Biopython versions. 

To use, your version of python must have the Biopython package installed. If using a version of python with Pip, simply run the following command in the command line:

```
pip install biopython
```
If using Jupiter notebooks, biopython should be pre-installed, but you can check by running 
```
conda install -c conda-forge biopython
```

For more information about Installation, see the repository [here](https://biopython.org/)

# Explanation of Functions
This repository contains various useful functions. I have included some documentation about each of the functions below:

## IO functions
Allows for simple opening and saving of SeqRecord files as standard .gbk files

```
open_SeqRecord_obj(path)
```
Open a .gbk file and return a SeqRecord object

```
save_SeqRecord_obj(SeqRecord_obj,
                   path)
```
Save a SeqRecord object to a .gbk file

## Base Functions
These are simple functions to search and edit SeqRecord objects that will be abstracted upon later

```
find_sequence(SeqRecord_obj,
              search_seq,
              forward_only=False,
              reverse_only=False):
```
Returns all locations of your exact search sequence. If directionality is important, use the toggles.

```
remove_sequence(SeqRecord_obj,
                start_input,
                end_input,
                maintain_ambiguious_features=True,
                supress_warnings=False):
```
Will remove a sequence based on a particular region within the sequence. This function takes a SeqRecord_obj in the form of a SeqRecord and will remove the sequence between the start and end input location. By default, it will not remove features that are edited, but can remove such features if maintain_ambiguious_features=False. It will warn users about editing such features unless supress_warnings=True. This is a base function that is abstracted upon in other functions.

**Note**: This function is strand agnostic. Therefore start must be before end

```
add_sequence(SeqRecord_obj,
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
```
Add a sequence to a specified start location within a SeqRecord_obj as a SeqRecord object. It can also make that newly added sequence a feature (as default). If a feature is made, passthrough arguments can be used to specify feature name, feature type, and feature color. the feature_qualifiers are a set of other data that could add extra information that different visualizers use to interpret the feature, but is unnecessary in most cases. The same consideration about editing sequence within other features is here (see above feature_warning function). This is a base function that is abstracted upon with other functions.

**Note**: this function is not strand agnostic. by default it adds your sequence in the forward direction, but if strand=-1 the sequence will be added in reverse at the same start location.

## Feature Functions
These functions allow you to manipulate features within a SeqRecord object. It is important to point out that features are an independent object, but SeqRecord objects often contain a list of features annotating the sequence. 
```
identify_feature(SeqRecord_obj,
                 feature):
```
finds a list of things that could be your feature. If you give it a string, it will find all features with that name. If you give it an individual feature, it will make sure that feature is in your SeqRecord_obj and return a list containing just that feature. If it is a list of feature objects already, it will just return that same list.
```
add_translation_to_feature(SeqRecord_obj,
                           Feature_obj,
                           supress_warnings=False):
```
Will add a translation to a feature object in a format understood by snapgene. Note this is not done in place, so make sure to add the feature back to your SeqRecord_obj
```
add_feature(SeqRecord_obj,
            start,
            end,
            feature_name='N/A',
            feature_type='CDS',
            feature_qualifiers_input=OrderedDict(),
            feature_color='#FF0000',
            strand=1,
            supress_warnings=False):
```
Adds a feature annotation to a SeqRecord_obj which is a SeqRecord object. This feature is added between the start and end location specified. Feature name, type, and color can all be specified. the feature_qualifiers can be used to pass an ordered dictionary containing extra information about the feature directly to the feature. This extra information could be used by other programs to glen more information about that feature, but is largely unnecessary. If the feature_type is a coding sequence ('CDS'), a translation is included. This is a base function that is abstracted upon in other functions here.

**Note**: this function does not alter sequence, it simply adds a feature annotation to existing sequence.
```
remove_feature(SeqRecord_obj,
               feature):
```
Removes a feature by name or by direct passing of the feature object.

## Abstracted Functions
These functions are just specific implementations of the functions listed above. These allow you to do some high-level functions like replace one sequence with another by combining both addition and removal of sequences

```
add_feature_from_sequence(SeqRecord_obj,
                          search_seq,
                          feature_name='N/A',
                          feature_type='CDS',
                          feature_qualifiers_input=OrderedDict(),
                          feature_color='#FF0000',
                          forward_only=False,
                          reverse_only=False,
                          supress_warnings=False):
```
adds feature using a search sequence. If multiple copies of the feature are present, it will add multiple copies. This is an abstraction of the add_feature function.
```
delete_seq_by_features(SeqRecord_obj,
                       feature_input,
                       maintain_ambiguious_features=True,
                       supress_warnings=False):
```
Remove a feature and the underlying sequence. This is an abstraction of the remove_sequence.
```
replace_sequence_by_locus(SeqRecord_obj,
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
```
Will replace a sequence specified by the location with a new input sequence This is an abstraction of the remove_sequence and add_sequence functions.
```
replace_sequence_by_feature(SeqRecord_obj,
                            target_feature,
                            sequence,
                            make_feature=False,
                            feature_name='N/A',
                            feature_type='CDS',
                            feature_qualifiers=OrderedDict(),
                            feature_color='#FF0000',
                            maintain_ambiguious_features=True,
                            supress_warnings=False):
```
Will take a target feature and replace all instances of that target feature with a new sequence. this is an abstraction of replace_sequence_by_locus.
```
translate_CDS(SeqRecord_obj,
              supress_warnings=False):
```
Will iterate through a SeqRecord_obj to add translations to all the features with type "CDS"
#Test
After import, run the following code in a folder with the test_in.gbk file included in this repository to ensure that everything is processing correctly 
```
test()
```
This test function ensures that all various functions work properly. it must be run in a folder with test_in.gbk file.

**Note**: The expected behavior includes warning about editing features in an ambiguous way. Their presence with the production of new gbk files is a successful run.
