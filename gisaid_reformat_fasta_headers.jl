#=
gisaid_reformat_fasta_headers.jl
    Reformat fasta headers to make gisaid_epi_isl record id.
    New fasta header contains country of exposure, lineaage, and 
    gisaid_epi_isl as elements 4 through 6 of description.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_04_05
=#

using Pkg
Pkg.activate(".")
Pkg.instantiate()
using ArgParse
using CSV
using DataFrames
using FASTX

function GetArgs()
    # from https://argparsejl.readthedocs.io/en/latest/argparse.html
    s = ArgParseSettings(description = "Reformat FASTA headers to include lineage, isl id, and country of exposure.",
                         version = "1.0",
                         add_version = true)

    @add_arg_table s begin
        "--input_file"
        required = true
        help = "Input file (required). Fasta file from GISAID."
        arg_type = String

        "--meta_file"
        required = true
        help = "Meta file (required). GISAID meta file corresponding to FASTA file."
        arg_type = String

        "--output_file"
        required = false
        help = "Output FASTA file. (default: write to screen)"
        arg_type = String

        "--replace"
        help = "Replace U's with T's in sequence."
        action = :store_true

        "--n"
        required = false
        default = 20
        help = "Minimum allowed number of N's in a row. (default = 20)"
        arg_type = Int
    end

    return parse_args(s)
end

function main()
    args = GetArgs()
    in_file = args["input_file"]
    meta_file = args["meta_file"]
    out_file = args["output_file"]
    replace_u = args["replace"]
    min_Ns = args["n"]

    ns = repeat("N", min_Ns)

    if ! isnothing(out_file)
        f = open(FASTA.Writer, out_file)
    else
        f = stdout
    end

    df = DataFrame(CSV.File(meta_file, delim = '\t'))
    
    reader = open(FASTA.Reader, in_file)
    for rec in reader
        seq = FASTA.sequence(rec)
        if ! isnothing(findfirst(ns, convert(String, seq)))
            println(stderr, "Too many N's: ", FASTA.identifier(rec))
            continue
        end

        # BioJulia handes description differntly from BioPython
        ident = FASTA.identifier(rec)
        desc = FASTA.description(rec)
        if ! isnothing(desc)
            ident = ident * " " * desc
        end

        items1 = split(ident, "|")
        items2 = split(items1[1], "/")
        new_id = join(items2[2:end], "/")
        new_id = replace(new_id, " " => "")
        new_id = replace(new_id, "'" => "-")
        new_id = replace(new_id, "," => "")

        df2 = df[df.strain .== new_id, :]
        if size(df2)[1] == 0
            println(stderr, "Missing ID in metafile: ", FASTA.identifier(rec))
            continue
        end

        new_desc = join([FASTA.identifier(rec)
                         df2[:, "pango_lineage"]
                         df2[:, "country_exposure"]
                         df2[:, "gisaid_epi_isl"]], "|")

        if ! replace_u
            new_seq = convert(String, seq)
        else
            new_seq = replace(convert(String, seq), "U" => "T")
        end

        new_rec = FASTA.Record(new_id, new_desc, new_seq)
        write(f, new_rec)
    end
    close(f)
end

main()
