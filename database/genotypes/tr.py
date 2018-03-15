from Bio.Seq import Seq
def u_to_t(recs):
    for rec in recs:
        rec.seq = rec.seq.back_transcribe()
        yield rec

