from pymol import cmd

def color_sel(color,selection,set_label_color=True):
    cmd.color(color,selection)
    if set_label_color:
        cmd.set("label_color", color, selection)

def renumber_residues(sel="polymer",start=1):
    """Renumber residues. Selection should be a single chain.
    Run from Pymol.
    See also this module for renumbering based on
    connectivity: https://pymolwiki.org/index.php/Renumber
    """
    class gen_resi_cl:
        def __init__(self,start):
            self.prev_new_resi = start
            self.prev_old_resi = None
        def __call__(self,cur_old_resi):
            if self.prev_old_resi is None:
                ret = self.prev_new_resi
            else:
                if self.prev_old_resi != cur_old_resi:
                    ret = self.prev_new_resi + 1
                else:
                    ret = self.prev_new_resi
            self.prev_old_resi = cur_old_resi
            self.prev_new_resi = ret
            print (ret,cur_old_resi)
            return ret

    space = dict(gen_resi=gen_resi_cl(start=start))
    cmd.alter(sel,"resi=gen_resi(resi)",space=space)
    cmd.rebuild()
    cmd.sort()

def remap_residues(index_file,sel="polymer"):
    import pandas as pd
    from string import ascii_uppercase as ins
    df = pd.read_csv(index_file,dialect="excel-tab")
    gr_df = df.groupby("id_from")
    for chain, group in gr_df:
        ind = {}
        prev_to_ind = None
        i_ins = 0
        for i_x,x in enumerate(group[["from_ind","to_ind"]].itertuples(index=False)):
            item = tuple(str(_) for _ in tuple(x))
            from_ind,to_ind = item
            ## runs of identical to_ind cause insertion codes
            if i_x > 0:
                if to_ind == prev_to_ind:
                    to_ind = to_ind+ins[i_ins]
                    i_ins += 1
                else:
                    i_ins = 0
            if from_ind not in ind:
                ind[from_ind] = to_ind
            prev_to_ind = item[1]

        space = dict(ind_map=ind)
        cmd.alter("{} and chain {}".format(sel,chain),
                  "resi=ind_map[resi]",space=space)
        cmd.rebuild()
        cmd.sort()

def set_labels():
    cmd.set("label_size",-4) #in Angstroms, scales with depth
    cmd.set("label_position",[0,0,10]) #[3,3,6] [0,0,5] bring forward
    cmd.set("label_outline_color","white")
    #cmd.set("label_shadow_mode", 2)
    cmd.set("label_font",14) #bolder font for scaling down bitmaps

def sel_str_resi(resi):
    return "resi "+"+".join([str(x) for x in resi])

def show_spheres(sel,sphere_transp=0,dot_radius=1,spikes=False,sphere_scale=1.):
    #cmd.set("dot_width",dot_width,"mut")
    cmd.set("sphere_transparency",sphere_transp,sel)
    cmd.set("sphere_scale", sphere_scale, sel)
    if spikes:
        ##in non-transparent protein, make mutations stand out
        ##by increasing sphere radius and adding spikes (dots)
        cmd.set("dot_density",1,sel)
        cmd.set("dot_radius",dot_radius,sel)
        cmd.show("dots",sel)
    cmd.show("spheres",sel)

def render():
    cmd.set("ray_opaque_background","off")
    cmd.set("antialias", 2)
    cmd.ray(2400,2400)


def paint_vars_f_Bin(vars_a,vars_b):
    cmd.set_view("""0.147443876, -0.624960124, -0.766605735, \
        0.453299403, -0.646183312, 0.613975227, \
        -0.879078567, -0.438027859, 0.188018516, \
        0.000000000, 0.000000000, -252.833984375, \
        30.680297852, 182.966445923, 30.680259705, \
        199.336242676, 306.331726074, -20.000000000""")
    cmd.show_as("cartoon","all")
    color_sel("grey","all")
    cmd.select("epi","resi 62-69+196-212")
    cmd.select("F2","resi 22-109")
    cmd.select("F1","resi 160-516")
    with open(vars_a,"r") as _:
        vars_a = _.readline().strip().split()
        freq_a = _.readline().strip().split()
    with open(vars_b,"r") as _:
        vars_b = _.readline().strip().split()
        freq_b = _.readline().strip().split()
    vars_sel = "resi "+"+".join(vars_a+vars_b)
    cmd.select("vars",vars_sel)
    cmd.show_as("spheres","vars")
    vars_a_sel = "resi "+"+".join(vars_a)
    vars_b_sel = "resi "+"+".join(vars_b)
    color_sel("orange","F1 and "+vars_a_sel)
    color_sel("green","F1 and "+vars_b_sel)
    color_sel("yellow","F2 and "+vars_a_sel)
    color_sel("cyan","F2 and "+vars_b_sel)
    color_sel("red","epi")
    set_labels()
    cmd.label("vars and name cb","'%s'%(resi,)")

def paint_vars_f(vars_file,show_var_labels=False):
    import pandas as pd
    ## Should labels be colored the same way as atoms?
    match_label_color = False
    cmd.bg_color("white")
    cmd.set("depth_cue",0)
    cmd.set("ray_trace_fog",0)
    cmd.show_as("cartoon","all")
    cmd.set("transparency", 0.5, "all")
    #cmd.show("surface", "not chain F")
    cmd.show("surface", "all")
    #chain_colors = ("gray40", "gray60", "gray80")
    chain_colors = ("palecyan","wheat","paleyellow")
    for chain,col in zip("DEF",chain_colors):
        sel = "chain {}".format(chain)
        color_sel(col,sel,set_label_color=match_label_color)
        cmd.set("cartoon_color",col,sel)
    cmd.select("epi","resi 62-69+196-212")
    cmd.select("F2","resi 22-109")
    cmd.select("F1","resi 160-516")
    cmd.select("near_epi","byres epi expand 10")
    color_sel("tv_green","epi",set_label_color=match_label_color) #"smudge"
    cmd.set("transparency", 0.3, "epi")

    vars_all = pd.read_csv(vars_file,dialect="excel-tab")
    for item in vars_all.itertuples():
        res_sel = "resi {}".format(item.pos)
        cmd.alter(res_sel,
                  "b={}".format(item.freq))
        rgb = (item.red,item.green,item.blue)
        rgb_str = "_".join([str(_) for _ in rgb])
        color_name = "var_freq_{}".format(rgb_str)
        cmd.set_color(color_name,rgb)
        color_sel(color_name,res_sel,set_label_color=match_label_color)
    vars_sel = "resi "+"+".join(vars_all.pos.astype(str))
    cmd.select("vars",vars_sel)
    #cmd.show_as("spheres","vars")
    #color_sel("red","vars",set_label_color=match_label_color)
    #cmd.spectrum("b","blue_red","vars",
    #             minimum=min(freq_all),
    #             maximum=max(freq_all),
    #             byres=1)
    if show_var_labels:
        set_labels()
        cmd.label("vars and name cb and near_epi and (chain F or epi)","'%s'%(resi,)")
    show_spheres("vars")
    show_spheres("vars and epi",spikes=True)
    ##TODO: see if I can use this script to create
    ##a legend for the gradient:
    ##https://pymolwiki.org/index.php/Spectrumbar

def run_paint_vars_f(show_var_labels=False):
    paint_vars_f(vars_file="RSV_F_vars.txt",show_var_labels=show_var_labels)

def make_pre_fusion_f():
    bu_pdb = "4jhw.pdb1.gz"
    id_pdb = bu_pdb.split(".")[0]
    fasta = id_pdb+"_f.fasta"
    pdb = id_pdb+"_f.pdb"
    chain_sel = "F"
    cmd.reinitialize()
    cmd.load(bu_pdb,"bu",2)
    for chain, state in (("F",2),("E",3),("D",4)):
        cmd.create(chain,"bu and chain {}".format(chain_sel),state,1)
        cmd.alter(chain,"chain='{}'".format(chain))
        cmd.rebuild()
        renumber_residues("{} and polymer and chain {}".format(chain,chain))
        cmd.rebuild()
    cmd.delete("bu")
    cmd.rebuild()
    cmd.save(fasta,"polymer",state=1,format="fasta")
    cmd.save(pdb,"polymer",state=1,format="pdb")

def make_post_fusion_f():
    bu_pdb = "3rrr.pdb1.gz"
    id_pdb = bu_pdb.split(".")[0]
    fasta = id_pdb+"_f.fasta"
    pdb = id_pdb+"_f.pdb"
    cmd.reinitialize()
    cmd.load(bu_pdb,"bu")
    for chain, chain1, chain2 in (("X","A","B"),
                                          ("Y","C","D"),
                                          ("Z","E","F")):
        cmd.alter("chain {}+{}".format(chain1,chain2),
                  "chain='{}'".format(chain))
        cmd.rebuild()
        renumber_residues("polymer and chain {}".format(chain))
        cmd.rebuild()
    for chain, chain0 in (("F","X"),("E","Y"),("D","Z")):
        cmd.alter("chain {}".format(chain0),"chain='{}'".format(chain))
        cmd.rebuild()
        cmd.create(chain,"bu and chain {}".format(chain))
        cmd.rebuild()
    cmd.delete("bu")
    cmd.rebuild()
    cmd.save(fasta,"polymer",format="fasta")
    cmd.save(pdb,"polymer",format="pdb")

def make_proteins_for_remapping():
    """Run from directory PDB with reference structures.
    This makes FASTA and reference PDB biounit structures
    numbered from 1.
    Run this first, then switch to R and generate remapping files
    with rsv_ali_to_index() and selection strings for the
    variations with rsv_pdb_selections().
    """
    make_pre_fusion_f()
    make_post_fusion_f()

def pnghack(filepath, width=1024, height=768):
    """Workaround if cmd.png() doesn't work"""
    cmd.viewport(width, height)  # Set resolution
    cmd.mpng(filepath, 1, 1)  # Use batch png mode with 1 frame only
    cmd.mplay()  # cmd.mpng needs the animation to 'run'

def remap_and_paint(pdb_root_sel=None,show_var_labels = False):
    """Run from pdb_out directory.
    Run this after running R code that creates PDB residue remapping
    files and Pymol selections files for the variants"""
    from os.path import join as pjoin
    from glob import glob
    pdb_dir = "../PDB"
    pdb_roots = ("3rrr_f","4jhw_f")
    index_root = "pdb_ref_F"
    subtypes = ("A","B")
    pdb_refs_all = {}
    for pdb_root in pdb_roots:
        out_dir=pdb_root
        pdb_inp = pjoin(pdb_dir,pdb_root+".pdb")
        pdb_refs = {}
        for subtype in subtypes:
            cmd.reinitialize()
            cmd.load(pdb_inp, "ref")
            index_file = glob(pjoin(out_dir,index_root+subtype+"_*.txt"))[0]
            remap_residues(index_file)
            pdb_ref = pjoin(out_dir,index_root+subtype+".pdb")
            cmd.save(pdb_ref,"polymer")
            pdb_refs[subtype] = pdb_ref
        with open(pdb_refs["A"]) as _a,\
            open(pdb_refs["B"]) as _b:
            assert _a.read() == _b.read(),\
                ("Expected remapped PDB files under {}"\
                    "to be identical for all subtypes").format(pdb_root)
        pdb_refs_all[pdb_root] = pdb_refs

    for subtype in subtypes:
        cmd.reinitialize()
        pdb_to = pdb_refs_all[pdb_roots[1]][subtype]
        pdb_from = pdb_refs_all[pdb_roots[0]][subtype]
        pdb_from_new = pdb_from.rsplit(".pdb",1)[0]+"_ali.pdb"
        cmd.load(pdb_from, "from")
        cmd.load(pdb_to, "to")
        cmd.align("from","to")
        cmd.save(pdb_from_new,"from")
        pdb_refs_all[pdb_roots[0]][subtype] = pdb_from_new

    pdb_roots_show = pdb_roots
    if not pdb_root_sel is None:
        pdb_roots_show = [ pdb_roots_show[_] for _ in pdb_root_sel ]

    for pdb_root in pdb_roots_show:
        out_dir = pdb_root
        pdb_refs = pdb_refs_all[pdb_root]
        cmd.reinitialize()
        cmd.load(pdb_refs["A"],"ref")
        run_paint_vars_f(show_var_labels=show_var_labels)
        cmd.select("none")
        cmd.set_view("""0.147443876, -0.624960124, -0.766605735, \
            0.453299403, -0.646183312, 0.613975227, \
            -0.879078567, -0.438027859, 0.188018516, \
            0.000000000, 0.000000000, -252.833984375, \
            30.680297852, 182.966445923, 30.680259705, \
            199.336242676, 306.331726074, -20.000000000""")
        cmd.refresh()
        return
        ray = 0
        cmd.png(pjoin(out_dir,"all_front.png"),width=1000,height=800,ray=ray)
        cmd.set_view("""0.632159054,    0.586480796,   -0.506372988,\
        0.166245610,   -0.740964770,   -0.650640249,\
        -0.756796241,    0.327127188,   -0.565905452,\
        0.000000000,    0.000000000, -265.902557373,\
        30.680297852,  182.966445923,   30.680259705,\
        171.127899170,  360.677093506,  -20.000000000""")
        cmd.refresh()
        cmd.png(pjoin(out_dir,"all_top.png"),width=1000,height=800,ray=ray)
