import itertools
import re

from Bio import SeqIO


def get_digests(sequence, nr_allowed_overdigestions):
    digests = re.findall(r".(?:(?<![KR](?!P)).)*", sequence)
    return "".join(list(itertools.chain(digests[0:nr_allowed_overdigestions + 1])))


def cleave_m(sequence):
    if sequence[0].upper() == "M":
        return sequence[1:]
    else:
        return sequence


def cleave_only(prot_db_name_in):
    # Parse the db instead of indexing it, because faster and we are going to append linearly anyway
    protein_db_in = SeqIO.parse(prot_db_name_in, format='fasta')

    # Remove the .fasta suffix and give the new db a proper name
    processed_comet_file_name = prot_db_name_in[:-6] + "_comet-cleaved-only.fasta"

    with open(processed_comet_file_name, "w+") as comet_processed_db_out:
        # Add the virtually cleaved proteins
        for original_seq_req in protein_db_in:
            cleaved_seq = cleave_m(str(original_seq_req.seq))
            comet_processed_db_out.write(
                ">%s %s\n%s\n" % (original_seq_req.id, original_seq_req.description, cleaved_seq))


def cleave_and_digest(prot_db_name_in, min_seq_len, nr_allowed_overdigestions):
    # Parse the db instead of indexing it, because faster and we are going to append linearly anyway
    protein_db_in = SeqIO.parse(prot_db_name_in, format='fasta')

    # Remove the .fasta suffix and give the new db a proper name
    processed_comet_file_name = prot_db_name_in[:-6] + "_comet-digested-cleaved.fasta"

    # # Make a copy of the original database
    # shutil.copy2(prot_db_name_in, processed_comet_file_name)

    # Append the virtually cleaved and trypsin digested proteins into the copied file
    with open(processed_comet_file_name, "a+") as comet_processed_db_out:

        # Add the virtually cleaved and trypsin digested proteins
        for original_seq_req in protein_db_in:
            processed_seq = get_digests(cleave_m(str(original_seq_req.seq)), nr_allowed_overdigestions)

            # Only save processed sequences with a minnimum length
            if len(str(processed_seq)) >= min_seq_len:
                # if nr_allowed_overdigestions == 0:
                #     custom_descr = "M-CLEAVED_TRYPSIN-DIGESTED_NO_OVERDIGESTION"
                # else:
                #     custom_descr = "M-CLEAVED_TRYPSIN-DIGESTED_%i_ALLOWED_OVERDIGESTIONS" % nr_allowed_overdigestions

                # Use .write() instead of SeqIO.write(), because faster
                # comet_processed_db_out.write(">%s %s\n%s\n" % ("", custom_descr, processed_seq))

                comet_processed_db_out.write(
                    ">%s %s\n%s\n" % (original_seq_req.id, original_seq_req.description, processed_seq))


if __name__ == '__main__':
    # cleave_and_digest("../" + "GitHub_test_files/uniprot-organism3A-mycobacterium-tuberculosis_proteins.fasta", 5, 1)
    # cleave_only("../" + "GitHub_test_files/Mycobacterium_marinum_proteome.fasta")
    cleave_only("../" + "GitHub_test_files/Mycobacterium_marinum_proteome_six_frame.fasta")
