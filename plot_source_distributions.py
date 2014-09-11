import matplotlib.pyplot as plt
from astropy.table import Table
from gammapy.catalog import FluxDistribution

tables = [('Population A', 'simulated_galaxy_1.fits'),
          ('Population B', 'simulated_galaxy_2.fits'),
	  ('1FHL Catalog', 'gll_psch_v07.fit')]

flux_distributions = []
for label, filename in tables:
    table = Table.read(filename)
    try:
	table.rename_column('Flux', 'flux')
    except:
	pass
    table.rename_column('flux', 'S')
    flux_distribution = FluxDistribution(table, label=label)
    flux_distributions.append(flux_distribution)

plt.figure(figsize=(4,3))
for flux_distribution in flux_distributions:
    flux_distribution.plot_integral_count()
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.ylabel(r'Integral Counts ($>$ Flux)', fontsize=10)
plt.xlabel(r'Flux/ph cm$^{-2}$ s$^{-1}$', fontsize=10)
#plt.ylim(ymin=0)
plt.legend(loc='lower left', fontsize='x-small')
plt.tight_layout()
plt.savefig('log_N_log_S.pdf')
