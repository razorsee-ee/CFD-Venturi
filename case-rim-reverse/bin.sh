# Helpers shared by the run scripts. They edit single values in place with sed,
# so system/ files stay identical to the git version after Allclean.

# setControl <keyword> <value>   e.g. setControl endTime 2500 ; setControl stopAt writeNow
setControl()
{
    sed -i "s/^\($1[[:space:]]\{1,\}\)[^;]*;/\1$2;/" system/controlDict
}

# residualControl on|off  (off: the solver ignores the convergence criteria)
residualControl()
{
    if [ "$1" = off ]
    then
        sed -i 's/^\([[:space:]]*\)residualControl$/\1residualControlOff/' system/fvSolution
    else
        sed -i 's/^\([[:space:]]*\)residualControlOff$/\1residualControl/' system/fvSolution
    fi
}

# latestTime: last saved iteration (decomposed or not)
latestTime()
{
    if [ -d processor0 ]
    then
        foamListTimes -processor -latestTime 2>/dev/null | tail -1
    else
        foamListTimes -latestTime 2>/dev/null | tail -1
    fi
}
