set "proteowizardversion=3.0.9844"
for %%F in (%~dp0\*.mzid) do (
"C:\Program Files\ProteoWizard\ProteoWizard %proteowizardversion%\idconvert.exe" %%F --pepXML -v -e .pep.xml -o %~dp0 >> %%F.log 2>&1
)
pause 

