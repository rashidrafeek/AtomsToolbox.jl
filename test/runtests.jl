using Test
using AtomsToolbox
using AtomsBase, Unitful
using AtomsIO
using StaticArrays: SVector

@testset "AtomsToolbox.jl" begin
    # Si primitive crystal 
    # Obtained using `find_primitive` in pymatgen of materials project mp-149
    frac_coords = [[0, 0, 0], [0.75, 0.75, 0.75]]
    lattice = [[3.33357328, 0.        , 1.92463943],
                [1.11119109, 3.14292303, 1.92463943],
                [0.        , 0.        , 3.84927886]] .* u"Å"
    symbols = fill(:Si, 2)
    si_primitive = periodic_system(symbols .=> frac_coords, lattice; fractional=true)
    
    si = load_system("data/Si.cif")
    
    @testset "getters" begin
        @testset "distance" begin
            @test distance(si_primitive, 1, 2) ≈ 2.35719227u"Å"
            @test distance(si_primitive, 1, 2; pbc=false) ≈ 7.07157681u"Å"

            box = cell_vectors(si_primitive)
            for i in 1:3
                @test distance(si_primitive, zeros(3)u"Å", box[i]) ≈ 0.0u"Å" atol=1e-8u"Å"
            end
        end

        @testset "angle" begin
            @test angle(si, 3, 1, 5) ≈ 60.0u"°"
            @test angle(si, 1, 6, 7) ≈ 109.47122063u"°"
            @test angle(si, 3, 2, 5) ≈ 109.47122063u"°"
            @test angle(si, 3, 2, 5; pbc=false) ≈ 58.51784589u"°"
        end

        @testset "dihedral" begin
            @test dihedral(si, 1,5,2,6; pbc=false) ≈ 0.0u"°" atol=1e-8u"°"
            @test dihedral(si, 5,6,1,4) ≈ -60.0u"°" atol=1e-8u"°"
            @test dihedral(si, 5,6,1,4; pbc=false) ≈ 30.0u"°" atol=1e-8u"°"
        end

        @testset "distance_matrix" begin
            @test isapprox(distance_matrix(si_primitive), [0.0 2.35719; 2.35719 0.0]u"Å"; atol=1e-5u"Å")
        end

        @testset "cell_parameters" begin
            @test all(cell_parameters(si_primitive) .≈ vcat(fill(3.84927886u"Å", 3), fill(60.0u"°", 3)))
            @test all(cell_parameters(si) .≈ vcat(fill(5.44370237u"Å", 3), fill(90.0u"°", 3)))
        end

        @testset "cell_lengths" begin
            @test all(cell_lengths(si_primitive) .≈ fill(3.84927886u"Å", 3))
            @test all(cell_lengths(si) .≈ fill(5.44370237u"Å", 3))
        end

        @testset "cell_angles" begin
            @test all(cell_angles(si_primitive) .≈ fill(60.0u"°", 3))
            @test all(cell_angles(si) .≈ fill(90.0u"°", 3))
        end

        @testset "volume" begin
            @test volume(si_primitive) ≈ 40.32952680u"Å^3"
            @test volume(si) ≈ 161.31810712u"Å^3"
        end

        @testset "scaled_position" begin
            @test scaled_position(si_primitive) ≈ [[0, 0, 0], [0.75, 0.75, 0.75]]
        end
    end

    @testset "PBC" begin
        # Define a simple cubic cell with side length 3 Å.
        box = (
            [3.0, 0.0, 0.0]u"Å", 
            [0.0, 3.0, 0.0]u"Å", 
            [0.0, 0.0, 3.0]u"Å"
        )

        # Create a small system with two atoms:
        #    - First at [0, 0, 0] Å
        #    - Second at [2.5, 2.5, 2.5] Å
        # We expect pbc_shortest_vector to produce [-0.5, -0.5, -0.5] Å 
        # as the minimal image distance from the first to second atom.
        atom1 = Atom(:X, [0.0, 0.0, 0.0]u"Å")
        atom2 = Atom(:X, [2.5, 2.5, 2.5]u"Å")

        # Build a FlexibleSystem (from AtomsBase) with periodic boundary conditions.
        test_sys = FlexibleSystem(
            [atom1, atom2]; 
            cell_vectors = box, 
            periodicity = (true, true, true)
        )

        @testset "Two-atom cubic system" begin
            vec = AtomsToolbox.pbc_shortest_vector(test_sys, 1, 2)
            @test vec ≈ SVector(-0.5, -0.5, -0.5)u"Å" atol=1e-10u"Å"

            # Also test the version that accepts explicit positions:
            pos1 = [0.0, 0.0, 0.0]u"Å"
            pos2 = [2.5, 2.5, 2.5]u"Å"
            vec2 = AtomsToolbox.pbc_shortest_vector(test_sys, pos1, pos2)
            @test vec2 ≈ SVector(-0.5, -0.5, -0.5)u"Å" atol=1e-10u"Å"
        end
    end
end
