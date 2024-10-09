### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ da3158ce-84dd-11ee-0e81-3516ae8782d9
begin
	# add and load general packages used in this notebook.
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor
end;

# ╔═╡ 1b5cb9ff-8487-4e75-b17b-3a373bb8a1b0
begin
	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end;

# ╔═╡ e3f3878b-67b8-4a28-b689-fac7a77d8136
md"""
## Part 0: Installing ACEpotentials
ACEpotentials normally requires Julia 1.9. For detailed installation instructions, see: https://acesuit.github.io/ACEpotentials.jl/dev/gettingstarted/installation/.

Warning: You may need to run the following command two or three times before it succeeds fully.
"""

# ╔═╡ 812e8964-f6a8-4574-9db1-a499f47e6aaf
md"""
##### ACEpotentials installation (requires Julia 1.9)
"""

# ╔═╡ 07540757-a430-456d-8053-db1f34842734
md"""
##### Checking the status of the installed packages!
"""

# ╔═╡ 14407001-9740-4f51-bcef-4047d4adcf34
Pkg.status()

# ╔═╡ 9af041ff-70f0-4893-b936-a7976a9c2c20
md"""
## Part 1: Basic dataset analysis
### Acessing Data
ACEpotentials provides quick access to several example datasets, which can be useful for testing. The following command lists these datasets.
"""

# ╔═╡ a2f772d5-7735-4b3d-a8ff-9f31b1eb3337
ACEpotentials.list_example_datasets()

# ╔═╡ 5701747b-761d-40fb-aa8a-fa5847b65731
md"""
Beginning by loading the tiny silicon dataset.
"""

# ╔═╡ e99b7521-99ff-42eb-a879-4418cd53806c
md"""
These data were taken from a larger set published with:

A. P. Bartók, J. Kermode, N. Bernstein, and G. Csányi, Machine Learning a General-Purpose Interatomic Potential for Silicon, Phys. Rev. X 8, 041048 (2018)
"""

# ╔═╡ 9a13a2f0-94dd-44d7-89c7-5080460ffa3c
#Tiny dataset
Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");

# ╔═╡ 264a590c-537e-4b83-b193-1dedbc16b518
md"""
To illustrate the procedure for loading extended xyz data from a file, we download the larger dataset and load it.
"""

# ╔═╡ bf1b81ac-97eb-46ee-abfb-f0b19cad2a3c
begin
	#Big dataset
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
end;

# ╔═╡ e4b8188f-48ae-47e9-af20-c870cddd588d
md"""
Assessing the dataset sizes.
"""

# ╔═╡ d61552f2-bf22-4fca-8622-8a001698aa41
begin
	println("The tiny dataset has ", length(Si_tiny_dataset), " structures.")
	println("The large dataset has ", length(Si_dataset), " structures.")
end

# ╔═╡ 63f15e27-f186-44f8-9812-6139240c21c6
md"""
Creating arrays containing the config_type for each structure in the datasets. Afterwards, we count the configurations of each type.
"""

# ╔═╡ 8ee22d81-db9e-4985-a0fe-bf2f4184e10d
begin
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	config_types = [at.data["config_type"].data for at in Si_dataset]
	
	function count_configs(config_types)
	    config_counts = [sum(config_types.==ct) for ct in unique(config_types)]
	    config_dict = Dict([ct=>cc for (ct,cc) in zip(unique(config_types), config_counts)])
	end;
end

# ╔═╡ 0093e684-dcba-43a1-bbf9-b43341b49adf
config_types 

# ╔═╡ 8a33d459-6fd1-44e4-986f-e51b4c852fba
begin
	#Counting Conf. of Tiny dataset
	println("There are ", length(unique(config_types_tiny)), " unique config_types "*
	        "in the tiny dataset:")
	display(count_configs(config_types_tiny))
end

# ╔═╡ cb3f1c8d-1042-43bf-a357-9562b140f743
begin
	#Counting Conf. of Big dataset
	println("There are ", length(unique(config_types)), " unique config_types "*
	        "in the full dataset:")
	display(count_configs(config_types))
end

# ╔═╡ 62d46ac2-e194-4730-8bbe-bc652f8056ee
md"""
### Analysis of Basis states
Two basic distributions which indicate how well the data fills space are the radial and angular distribution functions. We begin with the radial distribution function, plotting using the histogram function in Plots.jl. For the RDF we add some vertical lines to indicate the distances and first, second neighbours and so forth to confirm that the peaks are in the right place.
"""

# ╔═╡ 7890b7b8-2bb4-440e-8ff9-f01ba316dc0a
md"""
###### We can see the larger dataset clearly has a better-converged radial distribution function.
"""

# ╔═╡ 0b79e0b5-3a83-4431-ae9a-81950c457c60
begin
	r_cut = 6.0
	
	rdf_tiny = ACEpotentials.get_rdf(Si_tiny_dataset, r_cut; rescale = true)
	plt_rdf_1 = histogram(rdf_tiny[(:Si, :Si)], bins=150, label = "rdf",
	                      title="Si_tiny_dataset", titlefontsize=10,
	                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
	                      xlims=(1.5,6), size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:Si)*[1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black")
	
	rdf = ACEpotentials.get_rdf(Si_dataset, r_cut; rescale = true);
	plt_rdf_2 = histogram(rdf[(:Si, :Si)], bins=150, label = "rdf",
	                      title="Si_dataset", titlefontsize=10,
	                      xlabel = L"r[\AA]", ylabel = "RDF", yticks = [],
	                      xlims=(1.5,6), size=(400,200), left_margin = 2Plots.mm)
	vline!(rnn(:Si)*[1.0, 1.633, 1.915, 2.3, 2.5], label = "r1, r2, ...", lw=5, color = "black")
	
	plot(plt_rdf_1, plt_rdf_2, layout=(2,1), size=(600,400))
end

# ╔═╡ 486c7dfa-10d9-4587-931e-45f4556d4477
md"""
For the angular distribution function, we use a cutoff just above the nearest-neighbour distance so we can clearly see the equilibrium bond-angles. In this case, the vertical line indicates the equilibrium bond angle.
"""

# ╔═╡ f3cba55e-19ee-445f-8caf-167483f3bdc0
begin
	r_cut_adf = 1.25 * rnn(:Si)
	eq_angle = 1.91 # radians
	adf_tiny = ACEpotentials.get_adf(Si_tiny_dataset, r_cut_adf);
	plt_adf_1 = histogram(adf_tiny, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "Si_tiny_dataset", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "109.5˚", lw=5, color = "black")
	
	adf = ACEpotentials.get_adf(Si_dataset, r_cut_adf);
	plt_adf_2 = histogram(adf, bins=50, label = "adf", yticks = [], c = 3,
	                    title = "Si_dataset", titlefontsize = 10,
	                    xlabel = L"\theta", ylabel = "ADF",
	                    xlims = (0, π), size=(400,200), left_margin = 2Plots.mm)
	vline!([ eq_angle,], label = "109.5˚", lw=5, color = "black")
	
	plot(plt_adf_1, plt_adf_2, layout=(2,1), size=(600,400))
end

# ╔═╡ d060defc-3881-4d13-9f11-158bf103ef15
md"""
For later use, we define a function that extracts the energies stored in the silicon datasets.
"""

# ╔═╡ 5b4944fd-411a-45ec-a97b-8abc38febf50
begin
	function extract_energies(dataset)
	    energies = []
	    for atoms in dataset
	        for key in keys(atoms.data)
	            if lowercase(key) == "dft_energy"
	                push!(energies, atoms.data[key].data/length(atoms))
	            end
	        end
	    end
	    return energies
	end;
	
	Si_dataset_energies = extract_energies(Si_dataset)
	
	GC.gc()
end

# ╔═╡ 5ddb8936-8f3b-49d9-9231-bf13d2044551


# ╔═╡ Cell order:
# ╠═e3f3878b-67b8-4a28-b689-fac7a77d8136
# ╠═da3158ce-84dd-11ee-0e81-3516ae8782d9
# ╠═812e8964-f6a8-4574-9db1-a499f47e6aaf
# ╠═1b5cb9ff-8487-4e75-b17b-3a373bb8a1b0
# ╠═07540757-a430-456d-8053-db1f34842734
# ╠═14407001-9740-4f51-bcef-4047d4adcf34
# ╟─9af041ff-70f0-4893-b936-a7976a9c2c20
# ╠═a2f772d5-7735-4b3d-a8ff-9f31b1eb3337
# ╟─5701747b-761d-40fb-aa8a-fa5847b65731
# ╟─e99b7521-99ff-42eb-a879-4418cd53806c
# ╠═9a13a2f0-94dd-44d7-89c7-5080460ffa3c
# ╟─264a590c-537e-4b83-b193-1dedbc16b518
# ╠═bf1b81ac-97eb-46ee-abfb-f0b19cad2a3c
# ╟─e4b8188f-48ae-47e9-af20-c870cddd588d
# ╠═d61552f2-bf22-4fca-8622-8a001698aa41
# ╟─63f15e27-f186-44f8-9812-6139240c21c6
# ╠═8ee22d81-db9e-4985-a0fe-bf2f4184e10d
# ╠═0093e684-dcba-43a1-bbf9-b43341b49adf
# ╠═8a33d459-6fd1-44e4-986f-e51b4c852fba
# ╠═cb3f1c8d-1042-43bf-a357-9562b140f743
# ╠═62d46ac2-e194-4730-8bbe-bc652f8056ee
# ╠═7890b7b8-2bb4-440e-8ff9-f01ba316dc0a
# ╠═0b79e0b5-3a83-4431-ae9a-81950c457c60
# ╠═486c7dfa-10d9-4587-931e-45f4556d4477
# ╠═f3cba55e-19ee-445f-8caf-167483f3bdc0
# ╠═d060defc-3881-4d13-9f11-158bf103ef15
# ╠═5b4944fd-411a-45ec-a97b-8abc38febf50
# ╠═5ddb8936-8f3b-49d9-9231-bf13d2044551
