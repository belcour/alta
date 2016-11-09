!include EnvVarUpdate.nsh
!include "MUI2.nsh"
!include "winmessages.nsh"


# Name of the installer
Name "ALTA"
Outfile "ALTA.exe"

InstallDir "$DESKTOP\ALTA"

;!define ALTADIR "${__FILEDIR__}\..\.."
!define ALTADIR "..\.."

# Request application privileges for Windows Vista
RequestExecutionLevel user


# Pages
  !define MUI_WELCOMEPAGE_TITLE "Welcome to ALTA Installer"

  !define MUI_HEADERIMAGE
  !define MUI_HEADERIMAGE_BITMAP "alta-logo-installer.bmp"

  !define MUI_ABORTWARNING
  !define MUI_ABORTWARNING_CANCEL_DEFAULT

  !define MUI_ICON "alta-logo-installer.ico"
  #!define MUI_WELCOMEFINISHPAGE_BITMAP "alta-logo-installer-welcome.bmp"

  !insertmacro MUI_PAGE_WELCOME


  !insertmacro MUI_PAGE_LICENSE "..\..\LICENSE.txt"
#	!insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES




  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES

;Interface Settings

  !insertmacro MUI_LANGUAGE "English"


# Section start
Section "ALTA" SecMain

  SetOutPath $INSTDIR\bin
  File "${ALTADIR}\build\softs\*.exe"

  SetOutPath $INSTDIR\lib
  File "${ALTADIR}\build\core\core.lib"

  SetOutPath $INSTDIR\plugins
  File "${ALTADIR}\build\plugins\nonlinear*.dll"
  File "${ALTADIR}\build\plugins\rational*.dll"
  File "${ALTADIR}\build\plugins\data*.dll"

  ; Python package is not compiling yet
  ;SetOutPath $INSTDIR\python
  ;File "${ALTADIR}\build\python\alta.dll"

  # Update the ENVIROMNENT
  WriteRegStr HKCU "Environment" "ALTA_DIR"   '$INSTDIR'
  WriteRegStr HKCU "Environment" "ALTA_PLUGIN_PATH"   '$INSTDIR\plugins'
  ${EnvVarUpdate} $0 "PATH"     "A" "HKCU" '$INSTDIR\bin'
  ${EnvVarUpdate} $1 "PYTHONPATH" "A" "HKCU" '$INSTDIR\python'



  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"

  ; make sure windows knows about the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

#End Section
SectionEnd

;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;FILES

  Delete "$INSTDIR\Uninstall.exe"

  ;Remove the directory and all its content!
  RMDir /r "$INSTDIR"


  ;Remove ALTA_PLUGIN_PATH AND ALTA_DIR env. variable
  DeleteRegKey HKCU "ALTA_DIR"
  DeleteRegKey HKCU "ALTA_PLUGIN_PATH"


  ${un.EnvVarUpdate} $0 "PATH" "R" "HKCU" "$INSTDIR"
  ${un.EnvVarUpdate} $0 "PYTHONPATH" "R" "HKCU" "$INSTDIR\python"


  ; make sure windows knows about the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

SectionEnd
