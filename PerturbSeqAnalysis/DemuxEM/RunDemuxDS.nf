params.inDir //The input 10X outs directory
params.numDS=50000
params.inRNA="${params.inDir}/raw_feature_bc_matrix.h5" //The output of 10X, the raw h5 file including crispr info
params.demuxCode="$projectDir/bin/demuxEM" //The command for demuxEM
params.getCSVScript="$projectDir/scripts/GetCSV.DS.py" //python script to extract crispr info
params.cleanScript="$projectDir/scripts/CleanOutput.py" //python script to clean output
params.outdir="${params.inDir}/DemuxEM_DS" //output directory

workflow{
    crispCSV=GetCrisprCSV(params.getCSVScript,params.inRNA)
    demux=RunDemuxEM(params.inRNA,crispCSV,params.demuxCode)
    res=CleanDeMuxOutput(demux,params.cleanScript)
}

process GetCrisprCSV
{
    input:
    path 'GetCSV.DS.py'
    path 'raw_feature_bc_matrix.h5'

    output:
    path 'Crisp.csv'

    '''
    python GetCSV.DS.py raw_feature_bc_matrix.h5 Crisp.csv 50000
    '''
}

process RunDemuxEM
{
    input:
    path 'raw_feature_bc_matrix.h5'
    path 'Crisp.csv'
    path 'demuxEM'

    output:
    path 'output_demux.zarr.zip'

    '''
    demuxEM raw_feature_bc_matrix.h5 Crisp.csv output
    '''


}

process CleanDeMuxOutput
{

    publishDir "${params.outdir}", mode: 'copy'

    input: 
    path 'output_demux.zarr.zip'
    path 'CleanOutput.py'

    output:
    path 'Results.csv'

    '''
    python CleanOutput.py output_demux.zarr.zip Results.csv
    '''

}
