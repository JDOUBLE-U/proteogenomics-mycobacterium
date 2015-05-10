

# Keywords between dollar signs are substituted by subversion.
# The date and build will only tell the truth after a branch or tag,
# since different files in trunk will have been changed at different times
date ="$Date: 2012-07-02 20:04:37 -0700 (Mon, 02 Jul 2012) $".split()[1]
revision = "$Revision: 149 $".split()[1]


__version__ = '3.4' 


description = "CoreBio %s (%s)" % (__version__,  date)


