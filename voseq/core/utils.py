import itertools
import re

from django.conf import settings

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data.CodonTable import TranslationError
from degenerate_dna import Degenera

from stats.models import Stats
from . import exceptions


def get_voucher_codes(cleaned_data):
    """Processes list of voucher codes entered by users.

    It receives data from a **Form class** (`cleaned_data`) and makes sure that
    there are not duplicated voucher codes passed to the dataset builder.
    It also drops voucher codes if specified by users using the double dash:
    `--CP100-10`.

    Args:
        * `form.cleaned_data`: taken from a Form class

    Returns:
        tuple of voucher codes, no dupes, dropped unwanted.
    """
    voucher_codes = tuple()
    if cleaned_data['taxonset'] is not None:
        voucher_codes += tuple(cleaned_data['taxonset'].taxonset_list.splitlines())
    if cleaned_data['voucher_codes'] != '':
        voucher_codes += tuple(cleaned_data['voucher_codes'].splitlines())

    voucher_codes_clean = tuple()
    for i in voucher_codes:
        if re.search('^--', i):
            i_clean = re.sub('^--', '', i)
            voucher_codes_clean += (i_clean,)
        else:
            voucher_codes_clean += (i,)

    voucher_codes_set = tuple()
    for i in voucher_codes_clean:
        if i not in voucher_codes_set and i.strip() != '':
            voucher_codes_set += (i,)

    vouchers_to_drop = []
    for i in voucher_codes:
        if re.search('^--', i):
            vouchers_to_drop.append(re.sub('^--', '', i))

    voucher_codes_filtered = tuple()
    for i in voucher_codes_set:
        if i not in vouchers_to_drop:
            voucher_codes_filtered += (i,)
    return voucher_codes_filtered


def get_gene_codes(cleaned_data):
    """Processes list of gene codes entered by users.

    It receives data from a **Form class** (`cleaned_data`) and makes sure that
    there are not duplicated gene codes passed to the dataset builder.

    Args:
        * `form.cleaned_data`: taken from a Form class

    Returns:
        set of gene codes, sorted, no duplicates.
    """
    gene_codes = []
    if cleaned_data['geneset'] is not None:
        geneset_list = cleaned_data['geneset'].geneset_list.splitlines()
        for i in geneset_list:
            if i not in gene_codes:
                gene_codes.append(i)

    if cleaned_data['gene_codes']:
        for i in cleaned_data['gene_codes']:
            if i.gene_code not in gene_codes:
                gene_codes.append(i.gene_code)

    return tuple(sorted(gene_codes, key=str.lower))


def get_version_stats():
    """
    Returns version and database statistics for page footer.
    """
    version = settings.VERSION
    try:
        stats = Stats.objects.get(pk=1)
    except Stats.DoesNotExist:
        stats = ''

    return version, stats


def get_username(request):
    username = 'Guest'
    if request.user.is_authenticated():
        username = request.user.username
    return username


def clean_positions(a_list):
    if 'ALL' in a_list:
        return ['ALL']
    elif '1st' in a_list and '2nd' in a_list and '3rd' in a_list:
        return ['ALL']
    elif '1st' in a_list and '2nd' in a_list and '3rd' not in a_list:
        return ['1st-2nd']
    elif len(a_list) == 2:  # 1st and 3rd, or 2nd and 3rd
        raise exceptions.InadequateCodonPositions(
            "Cannot create dataset for only codon positions {0} and {1}.".format(a_list[0], a_list[1])
        )
    else:
        return a_list


def flatten_taxon_names_dict(dictionary):
    """Converts a dict to string suitable for FASTA object id

    Args:
        ``dictionary``: {'code': 'CP100-10', 'orden': 'Lepidoptera'. 'genus': 'Danaus'}

    Returns:
        Flattened as string: 'CP100-10_Lepidoptera_Danaus'

    """
    out = ''
    if 'code' in dictionary:
        out += dictionary['code'] + "_"

    if 'orden' in dictionary:
        out += dictionary['orden'] + "_"

    if 'superfamily' in dictionary:
        out += dictionary['superfamily'] + "_"

    if 'family' in dictionary:
        out += dictionary['family'] + "_"

    if 'subfamily' in dictionary:
        out += dictionary['subfamily'] + "_"

    if 'tribe' in dictionary:
        out += dictionary['tribe'] + "_"

    if 'subtribe' in dictionary:
        out += dictionary['subtribe'] + "_"

    if 'genus' in dictionary:
        out += dictionary['genus'] + "_"

    if 'species' in dictionary:
        out += dictionary['species'] + "_"

    if 'subspecies' in dictionary:
        out += dictionary['subspecies'] + "_"

    if 'author' in dictionary:
        out += dictionary['author'] + "_"

    if 'hostorg' in dictionary:
        out += dictionary['hostorg'] + "_"

    out_striped = re.sub('_+', '_', out)
    out_clean = re.sub('_$', '', out_striped)
    return out_clean.replace(" ", "_")


def chain_and_flatten(seqs):
    """Takes seq objects which only contain certain codon positions.

    Combines the two seq objects and returns another seq object.

    """
    out = []
    append = out.append

    my_chain = itertools.zip_longest(seqs[0], seqs[1])
    for i in itertools.chain.from_iterable(my_chain):
        if i is not None:
            append(i)

    return Seq(''.join(out))


def get_start_translation_index(gene_model, removed):
    start_translation = 0
    if int(gene_model['reading_frame']) == 1:
        if removed % 3 == 0:
            start_translation = 0
        if removed % 3 == 1:
            start_translation = 2
        if removed % 3 == 2:
            start_translation = 1
    elif int(gene_model['reading_frame']) == 2:
        if removed % 3 == 0:
            start_translation = 1
        if removed % 3 == 1:
            start_translation = 0
        if removed % 3 == 2:
            start_translation = 2
    elif int(gene_model['reading_frame']) == 3:
        if removed % 3 == 0:
            start_translation = 2
        if removed % 3 == 1:
            start_translation = 1
        if removed % 3 == 2:
            start_translation = 0
    else:
        raise exceptions.MissingReadingFrameForGene("Gene {0}".format(gene_model['gene_code']))
    return start_translation


def degenerate(gene_model, sequence, degen_translation):
    translation_start_position = get_start_translation_index(gene_model, removed=0) + 1
    bases_to_remove = translation_start_position - 1

    my_type = degen_translation
    if my_type.upper() == 'NORMAL':
        my_type = 'normal'

    dna = sequence[bases_to_remove:]
    res = Degenera(dna.upper(), gene_model['genetic_code'], my_type)
    res.degenerate()

    # put back the base that was excluded from the degeneration
    missing = sequence[:bases_to_remove].replace('?', 'N').upper()
    out = '{}{}'.format(missing, res.degenerated)

    return out


def gapped_translation(seq_obj, genetic_code):
    gap_indexes, sequence = get_gap_indexes(seq_obj)
    seq = Seq(sequence, generic_dna)

    ungapped_seq = seq.ungap('-')
    ungapped_seq = str(ungapped_seq).replace('?', 'N')
    ungapped_seq = Seq(ungapped_seq, generic_dna)

    translated_seq = ungapped_seq.translate(table=genetic_code)
    translated_seq_with_gaps = add_gaps_to_seq(translated_seq, gap_indexes)
    return str(translated_seq_with_gaps)


def add_gaps_to_seq(aa_sequence, gap_indexes):
    aa_seq_as_list = list(aa_sequence)

    number_of_question_marks_appended = 0
    for index in gap_indexes:
        new_index = index - number_of_question_marks_appended
        try:
            this_aa = aa_seq_as_list[new_index]
        except IndexError:
            this_aa = ''

        try:
            aa_seq_as_list[new_index] = '?' + this_aa
        except IndexError:
            aa_seq_as_list.append('?')
        number_of_question_marks_appended += 1
    return ''.join(aa_seq_as_list)


def get_gap_indexes(seq_obj):
    """If - is found not forming gap codons, it will be replaced by ? and
    the new sequence will be returned with this replacemen.
    """
    indexes_for_gaps_in_translated_sequence = []
    new_sequence = ''

    i = 0
    for index in range((len(seq_obj) // 3) + 1):
        j = i + 3
        tmp = str(seq_obj[i:j])
        if tmp.find('---') == 0:
            indexes_for_gaps_in_translated_sequence.append(index)
        elif '-' in tmp:
            tmp = tmp.replace('-', '?')
        new_sequence += tmp
        i += 3
    return indexes_for_gaps_in_translated_sequence, new_sequence


def strip_question_marks(seq):
    removed = 0
    seq = seq.upper()

    res = re.search('^(\?+)', seq)
    if res:
        removed += len(res.groups()[0])
    seq = re.sub('^\?+', '', seq)

    res = re.search('^(N+)', seq)
    if res:
        removed += len(res.groups()[0])
    seq = re.sub('^N+', '', seq)

    seq = re.sub('\?+$', '', seq)
    seq = re.sub('N+$', '', seq)
    return seq, removed
