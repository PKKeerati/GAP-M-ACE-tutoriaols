### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 15830db7-7ef6-4222-838d-1a86eff4ca85
begin
	using Pkg
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")
	Pkg.add("MultivariateStats")
	Pkg.add("Plots")
	Pkg.add("Suppressor")
	using LaTeXStrings, MultivariateStats, Plots, Printf, Statistics, Suppressor#

	Pkg.activate(".")
	Pkg.Registry.add("General")  # only needed when installing Julia for the first time
	Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry"))
	Pkg.add("ACEpotentials")
	using ACEpotentials
end

# ╔═╡ 9cc36068-1bbd-485b-a8d0-3837c54c7fa7
md"""
### Installing ACEpotentials
"""

# ╔═╡ 92b9901f-e3f6-4f89-b819-127e5f6b83f7
md"""
### Downloading dataset
"""

# ╔═╡ b4423e07-1b6e-4a77-ae50-d1f2ecd70d1a
begin
	Si_tiny_dataset, _, _ = ACEpotentials.example_dataset("Si_tiny");
	config_types_tiny = [at.data["config_type"].data for at in Si_tiny_dataset]
	
	
	download("https://www.dropbox.com/scl/fi/mzd7zcb1x1l4rw5eswxcd/gp_iter6_sparse9k.xml.xyz?rlkey=o4avtpkka6jnqn7qg375vg7z0&dl=0",
	         "Si_dataset.xyz");
	
	Si_dataset = read_extxyz("Si_dataset.xyz");
	config_types = [at.data["config_type"].data for at in Si_dataset]
end;

# ╔═╡ d6c59890-d1e2-4bb2-b281-db1ca4b667e0
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

# ╔═╡ ccab2500-d930-41ba-bfea-86ba9bf7026f
md"""
## Part 4: Committee models
ACEpotentials.jl can produce committee models using Bayesian linear regression. Such committees provide uncertainty estimates useful for active learning.

Recall our two silicon datasets. We begin by training a (relatively small) model on the tiny version.

Note the use of the BLR solver with a nonzero committee size.
"""

# ╔═╡ dc0eccf0-8541-11ee-2421-0713757269ef
begin
	model = acemodel(elements = [:Si,],
	                 order = 3,
	                 totaldegree = 12,
	                 Eref = [:Si => -158.54496821]);
	
	acefit!(model, Si_tiny_dataset;
	        solver = ACEfit.BLR(committee_size=50, factorization=:svd),
	        energy_key = "dft_energy", force_key = "dft_force",
	        verbose = false);
end

# ╔═╡ 6ef7f775-44bd-4ce2-be9d-29b0ea815448
md"""
Next we define a function which assesses model performance on the full silicon dataset.
"""

# ╔═╡ c2481c27-3bf4-474e-b0fa-4932ddcea13b
function assess_model(model, train_dataset)

    plot([-164,-158], [-164,-158]; lc=:black, label="")

    model_energies = []
    model_std = []
    for atoms in Si_dataset
        ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ene/length(atoms))
        push!(model_std, std(co_ene/length(atoms)))
    end
	
    rmse = sqrt(sum((model_energies-Si_dataset_energies).^2)/length(Si_dataset))
    mae = sum(abs.(model_energies-Si_dataset_energies))/length(Si_dataset)


    scatter!(Si_dataset_energies, model_energies;
             label="full dataset",
             title = @sprintf("Structures Used In Training:  %i out of %i\n", length(train_dataset), length(Si_dataset)) *
                     @sprintf("RMSE (MAE) For Entire Dataset:  %.0f (%.0f) meV/atom", 1000*rmse, 1000*mae),
             titlefontsize = 8,
             yerror = model_std,
             xlabel="Energy [eV/atom]", xlims=(-164,-158),
             ylabel="Model Energy [eV/atom]", ylims=(-164,-158),
             aspect_ratio = :equal, color=1)

    model_energies = [energy(model.potential,atoms)/length(atoms) for atoms in train_dataset]
    scatter!(extract_energies(train_dataset), model_energies;
             label="training set", color=2)

end;

# ╔═╡ 4609fcce-d2e8-4106-8865-56dc20f7b6a9
md"""
Applying this function to our current model yields
"""

# ╔═╡ 4c8e1eb4-aedc-4344-808f-5cfb6e37d607
assess_model(model, Si_tiny_dataset)

# ╔═╡ 200598e9-1d14-41b6-bf6b-0a3d9c7299f0
md"""
Clearly there is room to improve: the model-derived RMSE is 280 meV/atom for the full dataset. Moreover, the error bars show the standard deviation of the energies predicted by the commmittee, which are quite high for some data.

Next, we will define a function that augments the tiny dataset by adding structures for which the model is least confident.
"""

# ╔═╡ adf9c84d-21df-44c1-bb57-f80e8d09fb72
function augment(old_dataset, old_model; num=5)

    new_dataset = deepcopy(old_dataset)
    new_model = deepcopy(old_model)

    model_std = []
    for atoms in Si_dataset
        ene, co_ene = ACE1.co_energy(new_model.potential, atoms)
        push!(model_std, std(co_ene/length(atoms)))
    end
	#Choosing training points 
    for atoms in Si_dataset[sortperm(model_std; rev=true)[1:num]]
        push!(new_dataset, atoms)
    end
    @suppress acefit!(new_model, new_dataset;
            solver = ACEfit.BLR(committee_size=50, factorization=:svd),
            energy_key = "dft_energy", force_key = "dft_force",
            verbose = false);

    return new_dataset, new_model
end;

# ╔═╡ 1472e38e-5d35-4038-b07d-8d2dcab22c9b
md"""
The following applies this strategy, adding the five structures with the highest committee deviation.
"""

# ╔═╡ 7c387fde-81df-415e-863d-ec585d328900
begin
	new_dataset, new_model = augment(Si_tiny_dataset, model; num=5);
	assess_model(new_model, new_dataset)
end

# ╔═╡ 80f17d73-81e0-44c7-8c6b-005474116ba9
md"""
Already, there is notable improvement. The overall errors have dropped, and the predictions for the worst-performing structures are much improved.

Next, we perform four additional augmentation steps, adding twenty structures in total.
"""

# ╔═╡ 6a99ca1d-6f23-4ba5-8dc3-f23601b4b67b
begin
	for i in 1:4
	    @show i
	    new_dataset, new_model = augment(new_dataset, new_model; num=5);
	end
	assess_model(new_model, new_dataset)
end

# ╔═╡ 43440c1c-d9da-4278-aa83-737939597ccd
begin
	for i in 1:4
	    @show i
	    new_dataset, new_model = augment(new_dataset, new_model; num=10);
	end
	assess_model(new_model, new_dataset)
end

# ╔═╡ 56dbb733-9b84-47a5-9aaf-7cb6625ec485
md"""
Remarkably, although we are using only a small fraction (~3%) of the full dataset, our model now performs reasonably well.

Further iterations may improve on this result; however, a larger model is necessary to obtain extremely low errors.

Important: While this dataset filtering can be useful, the connection with active learning is crucial. Recall that we did not use the reference energies when selecting structures, only the committee deviation.
"""

# ╔═╡ dcd08f1d-04e2-4edb-a994-3d67164ee6de
md"""
## Model energy and STD
"""

# ╔═╡ 6c617be8-550e-43e7-8667-98878f8d96ce
begin
	model_energies = []
    model_std = []
    for atoms in Si_dataset
        ene, co_ene = ACE1.co_energy(model.potential, atoms)
        push!(model_energies, ene/length(atoms))
        push!(model_std, std(co_ene/length(atoms)))
    end
end

# ╔═╡ 5774584f-e7d3-4049-98d5-be5d0f7a6970
model_energies

# ╔═╡ 660ed92d-a653-47dd-a344-c6e5f7ea6163
begin
	new_model_energies = []
	new_model_std = []
	    for atoms in Si_dataset
	        ene, co_ene = ACE1.co_energy(new_model.potential, atoms)
	        push!(new_model_energies, ene/length(atoms))
	        push!(new_model_std, std(co_ene/length(atoms)))
	    end
end

# ╔═╡ a5e1caa4-5de6-48e2-ab4e-249362b67af6
new_model_energies

# ╔═╡ 725139f6-1552-4880-b17b-2861e5af4371
new_model_std

# ╔═╡ b25a235d-11f6-4dc5-8358-ccb356643890


# ╔═╡ 81992eb6-8fa0-473d-aa34-d8b54045b17d
GC.gc()

# ╔═╡ Cell order:
# ╟─9cc36068-1bbd-485b-a8d0-3837c54c7fa7
# ╠═15830db7-7ef6-4222-838d-1a86eff4ca85
# ╟─92b9901f-e3f6-4f89-b819-127e5f6b83f7
# ╠═b4423e07-1b6e-4a77-ae50-d1f2ecd70d1a
# ╠═d6c59890-d1e2-4bb2-b281-db1ca4b667e0
# ╟─ccab2500-d930-41ba-bfea-86ba9bf7026f
# ╠═dc0eccf0-8541-11ee-2421-0713757269ef
# ╟─6ef7f775-44bd-4ce2-be9d-29b0ea815448
# ╠═c2481c27-3bf4-474e-b0fa-4932ddcea13b
# ╟─4609fcce-d2e8-4106-8865-56dc20f7b6a9
# ╠═4c8e1eb4-aedc-4344-808f-5cfb6e37d607
# ╟─200598e9-1d14-41b6-bf6b-0a3d9c7299f0
# ╠═adf9c84d-21df-44c1-bb57-f80e8d09fb72
# ╟─1472e38e-5d35-4038-b07d-8d2dcab22c9b
# ╠═7c387fde-81df-415e-863d-ec585d328900
# ╟─80f17d73-81e0-44c7-8c6b-005474116ba9
# ╠═6a99ca1d-6f23-4ba5-8dc3-f23601b4b67b
# ╠═43440c1c-d9da-4278-aa83-737939597ccd
# ╟─56dbb733-9b84-47a5-9aaf-7cb6625ec485
# ╟─dcd08f1d-04e2-4edb-a994-3d67164ee6de
# ╠═6c617be8-550e-43e7-8667-98878f8d96ce
# ╠═5774584f-e7d3-4049-98d5-be5d0f7a6970
# ╠═660ed92d-a653-47dd-a344-c6e5f7ea6163
# ╠═a5e1caa4-5de6-48e2-ab4e-249362b67af6
# ╠═725139f6-1552-4880-b17b-2861e5af4371
# ╠═b25a235d-11f6-4dc5-8358-ccb356643890
# ╠═81992eb6-8fa0-473d-aa34-d8b54045b17d
