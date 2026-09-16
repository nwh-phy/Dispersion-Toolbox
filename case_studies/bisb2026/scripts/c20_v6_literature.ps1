param([Parameter(Mandatory=$true)][string]$OutputDirectory)
$dois=@('10.1038/s41586-025-08711-x','10.1021/acs.nanolett.5c03319','10.1038/s41699-017-0003-9','10.1103/PhysRevB.95.094301','10.1103/PhysRevB.101.205412','10.1021/acsnano.5c07482')
$records=@()
foreach($doi in $dois){
 $r=Invoke-RestMethod -Uri ('https://api.crossref.org/works/'+$doi)
 $m=$r.message
 $records += [pscustomobject]@{doi=$m.DOI;title=$m.title;authors=$m.author;published=$m.published;journal=$m.'container-title';volume=$m.volume;page=$m.page;source=('https://api.crossref.org/works/'+$doi);checked='2026-09-12'}
}
$records | ConvertTo-Json -Depth 10 | Set-Content -LiteralPath (Join-Path $OutputDirectory 'literature_metadata.json') -Encoding UTF8
