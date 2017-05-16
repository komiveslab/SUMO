for %%F in (%~dp0\*.raw) do (
C:\Inetpub\tpp-bin\msconvert %%F -v -z --mzXML -o %~dp0 >> %%F.log 2>&1
)
pause 