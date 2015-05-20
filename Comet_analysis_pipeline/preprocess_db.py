import itertools
import os
import re

from Bio import SeqIO


def get_digests(sequence, nr_allowed_overdigestions=2):
    # Default allowed overdigestion of Comet is 2 as well
    digests = re.findall(r".(?:(?<![KR](?!P)).)*", sequence)
    return "".join(list(itertools.chain(digests[0:nr_allowed_overdigestions + 1])))


def cleave_m(sequence):
    # Check may not be necessary
    if sequence[0].upper() == "M":
        return sequence[1:]
    else:
        return sequence


def cleave_m_only(prot_db_name_in, min_seq_len):

    print("Cleaving the M's of the protein database")

    # Parse the db instead of indexing it, because this is faster and we are going to append chronologically anyway
    protein_db_in = SeqIO.parse(prot_db_name_in, format='fasta')

    # Remove the .fasta suffix and give the new db an absolute path and a proper name
    procesed_db = os.getcwd() + "\\preprocessed_db_out\\" + prot_db_name_in[prot_db_name_in.rfind("/") + 1:-(
        len(".fasta"))] + "_comet_M_cleaved_only.fasta"

    with open(procesed_db, "w+") as comet_processed_db_out:
        # Add the virtually cleaved proteins
        for original_seq_req in protein_db_in:
            cleaved_seq = cleave_m(str(original_seq_req.seq))

            # Only save processed protein sequences with a given minnimum length
            if len(cleaved_seq) >= min_seq_len:
                comet_processed_db_out.write(
                    ">%s %s\n%s\n" % (original_seq_req.id, original_seq_req.description, cleaved_seq))

    return procesed_db


def cleave_m_and_digest(prot_db_name_in, min_seq_len):
    # Parse the db instead of indexing it, because this is faster and we are going to append chronologically anyway
    protein_db_in = SeqIO.parse(prot_db_name_in, format='fasta')

    # Remove the .fasta suffix and give the new db an absolute path and a proper name
    procesed_db = os.getcwd() + "\\preprocessed_db_out\\" + prot_db_name_in[prot_db_name_in.rfind("/") + 1:-(
        len(".fasta"))] + "_comet_digested_and_cleaved.fasta"

    # Append the virtually cleaved and trypsin digested proteins into the copied file
    with open(procesed_db, "a+") as comet_processed_db_out:

        # Add the virtually cleaved and trypsin digested proteins
        for original_seq_req in protein_db_in:
            processed_seq = get_digests(cleave_m(str(original_seq_req.seq)))

            # Only save processed protein sequences with a given minnimum length
            if len(str(processed_seq)) >= min_seq_len:
                # if nr_allowed_overdigestions == 0:
                #     custom_descr = "M-CLEAVED_TRYPSIN-DIGESTED_NO_OVERDIGESTION"
                # else:
                #     custom_descr = "M-CLEAVED_TRYPSIN-DIGESTED_%i_ALLOWED_OVERDIGESTIONS" % nr_allowed_overdigestions

                # Use .write() instead of SeqIO.write(), because faster
                # comet_processed_db_out.write(">%s %s\n%s\n" % ("", custom_descr, processed_seq))

                comet_processed_db_out.write(
                    ">%s %s\n%s\n" % (original_seq_req.id, original_seq_req.description, processed_seq))

    return procesed_db
