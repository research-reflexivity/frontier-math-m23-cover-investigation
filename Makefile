
.PHONY: paper verify verify-all verify-optimal verify-bridge verify-fano \
	hurwitz-pilot hurwitz-numerical hurwitz-multicentre hurwitz-acb \
	acb-fft hurwitz-acb-model verify-hurwitz-candidate \
	verify-hurwitz-models-candidate verify-hurwitz-marked-candidate \
	verify-hurwitz-branch-candidate \
	verify-hurwitz-maps-candidate verify-hurwitz-third-fiber-record \
	verify-hurwitz-third-fiber-exact reconstruct-hurwitz-maps \
	verify-hurwitz-monodromy-eliminant verify-hurwitz-branch-cycles \
	verify-hurwitz-relative-transporter verify-hurwitz-galois-closure \
	verify-hurwitz-local-23 verify-hurwitz-connector \
	verify-hurwitz-connector-a6 \
	certify-hurwitz-branch-cycles certify-degree-one-branch-cycles \
	hurwitz-monodromy-resultant \
	verify-hurwitz-tail-record verify-hurwitz-tail-geometry \
	hurwitz-tail-stability hurwitz-tail-model \
	reconstruct-canonical-quadric verify-magma-record verify-magma clean export-public

SAGE ?= sage
SINGULAR ?= Singular
MAGMA ?= magma
DOT_SAGE ?= /private/tmp/m23-cover-investigation-sage
HURWITZ_CLASS ?= 6
HURWITZ_TARGET ?= 2
HURWITZ_PRECISION ?= 256
GP_CERT = python3 scripts/run_gp_certificate.py

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	mkdir -p output/pdf
	cp paper/main.pdf output/pdf/m23-cover-investigation.pdf

verify: verify-all

hurwitz-pilot:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python scripts/compute_hurwitz_covers.py \
		--class-id 6 --terms 30 --samples 92

hurwitz-numerical:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python scripts/compute_hurwitz_covers.py \
		--all-classes --canonical --terms 120 --samples 299

hurwitz-multicentre:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python scripts/hurwitz_high_precision.py \
		--class-id 6 --terms 30 --samples 96 --canonical

hurwitz-acb:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python scripts/certify_hurwitz_acb.py \
		--class-id 6 --terms 30 --samples 96

acb-fft:
	cd scripts && DOT_SAGE=$(DOT_SAGE) $(SAGE) -python setup_acb_fft.py \
		build_ext --inplace

hurwitz-acb-model: acb-fft
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python scripts/compute_hurwitz_acb_model.py \
		--class-id 6 --terms 120 --samples 320 --precision 384 \
		--refine-rounds 2 --neumann-iterations 200

verify-hurwitz-tail-record:
	python3 verification/verify_hurwitz_tail_summary.py

verify-hurwitz-tail-geometry:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_tail_geometry.py

hurwitz-tail-stability:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		scripts/analyze_hurwitz_tail_stability.py \
		--class-id $(HURWITZ_CLASS) --low-terms 60 --samples 1280 \
		--precision 192 --certify-left-inverse

hurwitz-tail-model: acb-fft
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		scripts/compute_hurwitz_acb_model.py \
		--class-id $(HURWITZ_CLASS) --terms 700 --samples 1280 \
		--precision 1024 --refine-rounds 5 --neumann-iterations 400

verify-hurwitz-candidate:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_algebra_candidate.py

verify-hurwitz-models-candidate:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_canonical_models_candidate.py

verify-hurwitz-marked-candidate:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_marked_points_candidate.py

verify-hurwitz-branch-candidate:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_degree23_branch_candidate.py

verify-hurwitz-maps-candidate:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_degree23_maps_candidate.py

verify-hurwitz-monodromy-eliminant:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_monodromy_eliminant_candidate.py

verify-hurwitz-branch-cycles:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_branch_cycle_summary.py

verify-hurwitz-relative-transporter:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		verification/certify_hurwitz_relative_transporter.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_relative_transporter.py

verify-hurwitz-galois-closure:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_galois_closure.py
	gp -q verification/verify_hurwitz_galois_closure.gp

verify-hurwitz-local-23:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_local_23.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		verification/verify_hurwitz_pointed_23.py

verify-hurwitz-connector:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		notes/certify_ade_gluing_marker.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		notes/certify_wild_parameter_orientation.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		notes/certify_e8_tail_isomorphism.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		notes/certify_p23_special_deformation_datum.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/audit_lifted_trace_fano_affine_incidence.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/audit_fano_affine_odd_fixed_point_lemma.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/audit_raw_lifted_trace_node_boundary.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/audit_returned_normalizer_trace_parity.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/audit_fano_affine_incidence_bridge.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/explore_untagged_lifted_trace_pairing.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q \
		notes/explore_tagged_tame_boundary.g
	python3 notes/certify_pinched_tag_finite_identities.py

verify-hurwitz-connector-a6:
	M23_ADE_ONLY_A6=1 DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		notes/explore_p23_ade_deformation.py

certify-hurwitz-branch-cycles:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		scripts/certify_hurwitz_branch_cycles.py \
		--class-id $(HURWITZ_CLASS) --precision $(HURWITZ_PRECISION)

certify-degree-one-branch-cycles:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		scripts/certify_degree_one_branch_cycles.py \
		--precision $(HURWITZ_PRECISION)

hurwitz-monodromy-resultant:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		scripts/compute_hurwitz_monodromy_resultants.py \
		--target $(HURWITZ_TARGET) --output-dir /private/tmp/m23-hurwitz-monodromy

verify-hurwitz-third-fiber-record:
	python3 verification/verify_hurwitz_degree23_third_fiber_summary.py

verify-hurwitz-third-fiber-exact:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) \
		verification/verify_hurwitz_degree23_third_fiber.sage

reconstruct-hurwitz-maps:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python \
		scripts/reconstruct_hurwitz_degree23_maps.py --component all

verify-all: verify-optimal verify-bridge verify-fano verify-hurwitz-candidate \
	verify-hurwitz-models-candidate verify-hurwitz-marked-candidate \
	verify-hurwitz-branch-candidate \
	verify-hurwitz-maps-candidate verify-hurwitz-third-fiber-record \
	verify-hurwitz-monodromy-eliminant verify-hurwitz-branch-cycles \
	verify-hurwitz-relative-transporter verify-hurwitz-galois-closure \
	verify-hurwitz-local-23 verify-hurwitz-connector \
	verify-hurwitz-tail-record verify-hurwitz-tail-geometry verify-magma-record

verify-optimal:
	python3 scripts/render_fint_gp.py --check
	cd verification && gp -q verify_boundary_valuations.gp
	cd verification && SINGULAR="$(SINGULAR)" python3 verify_optimal_23_4.py
	cd verification && $(SINGULAR) --cpus=1 --threads=1 -q verify_nodes_mod31.sing
	python3 verification/verify_adjoint_mod31.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python verification/verify_sage.py
	python3 verification/verify_specialization_t2.py
	cd verification && gp -q verify_specialization_t2.gp
	python3 verification/verify_progression_319.py
	cd verification && gp -q verify_progression_319.gp
	python3 verification/verify_ramified_t3830.py
	cd verification && gp -q -s 4000000000 verify_ramified_t3830.gp

verify-bridge:
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python verification/verify_special_fibre_bridge.py

results:
	mkdir -p results

verify-fano: | results
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -python verification/verify_canonical_quadric.py
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_branch_cycle_description.g
	$(GP_CERT) verification/fano/certify_branch_cycle_description.gp
	$(GP_CERT) verification/fano/certify_special_fibre_fano_arithmetic.gp
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_special_fibre_fano_geometry.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_fano_flag_descent.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_local_flag_places.g
	$(GP_CERT) verification/fano/certify_local_flag_places.gp
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_special_fibre_field_tower.g
	$(GP_CERT) verification/fano/certify_special_fibre_field_tower.gp
	$(GP_CERT) verification/fano/verify_octic_splitting_degree.gp
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_vanishing_cycle_class.g
	$(GP_CERT) verification/fano/certify_vanishing_cycle_class.gp
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_special_fibre_ramification.g
	$(GP_CERT) verification/fano/certify_special_fibre_ramification.gp
	$(GP_CERT) verification/fano/certify_central_quadratic_ramification.gp
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_real_frame_cocycle.g
	$(GP_CERT) verification/fano/certify_real_frame_cocycle.gp
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_arithmetic_reflection_lift.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/certify_geometric_reflection_obstruction.g
	DOT_SAGE=$(DOT_SAGE) $(SAGE) -gap -A -q verification/fano/certify_gauss_prolongation_obstruction.g

reconstruct-canonical-quadric:
	python3 scripts/reconstruct_canonical_quadric_Q.py

verify-magma-record:
	python3 scripts/emit_magma_certificate.py --check
	python3 scripts/emit_canonical_quadric_magma_certificate.py --check
	python3 scripts/emit_hurwitz_degree23_magma_certificate.py --check
	python3 scripts/emit_hurwitz_degree23_magma_geometry_certificate.py --check
	python3 verification/verify_magma_summary.py
	python3 verification/verify_geometric_reflection_obstruction_magma_summary.py
	python3 verification/verify_hurwitz_degree23_branch_magma_summary.py
	python3 verification/verify_hurwitz_degree23_geometry_magma_summary.py
	python3 verification/verify_hurwitz_galois_closure_magma_summary.py
	python3 verification/verify_hurwitz_local_23_magma_summary.py

verify-magma: verify-magma-record
	$(MAGMA) -b verification/verify_optimal_23_4.m
	$(MAGMA) -b verification/certify_gauss_prolongation_obstruction.m
	$(MAGMA) -b verification/certify_geometric_reflection_obstruction.m
	$(MAGMA) -b verification/verify_canonical_quadric.m
	$(MAGMA) -b verification/verify_hurwitz_degree23_branch.m
	$(MAGMA) -b verification/verify_hurwitz_degree23_geometry.m
	$(MAGMA) -b verification/verify_hurwitz_galois_closure.m
	$(MAGMA) -b verification/verify_hurwitz_local_23.m
	$(MAGMA) -b notes/certify_fano_affine_odd_fixed_point_lemma.m
	$(MAGMA) -b notes/certify_pinched_tag_finite_identities.m
	$(MAGMA) -b notes/certify_wild_parameter_orientation.m

export-public:
	@test -n "$(EXPORT_DIR)" || (echo "Set EXPORT_DIR to the independent public repository."; exit 2)
	python3 scripts/export_public.py "$(EXPORT_DIR)"

clean:
	cd paper && latexmk -C main.tex
