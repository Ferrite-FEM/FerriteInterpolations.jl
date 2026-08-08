# L2 Piola mapping for scalar interpolations.
#
# Ferrite ships identity, covariant Piola and contravariant Piola mappings;
# DefElement additionally uses the "L2 Piola" mapping (e.g. for dPc), which
# pushes scalar reference values forward as
#
#     N(x) = N̂(ξ) / det J.
#
# This is the mapping that makes the L2 end of a de Rham complex commute on
# non-affinely mapped cells: div(contravariant-Piola v) = (1/det J) d̂iv(v̂),
# so the divergence of an H(div) element lands exactly in the L2-Piola-mapped
# scalar space. On affine cells det J is constant and the spanned space
# coincides with the identity-mapped one.
#
# Ferrite's FEValues machinery is extensible by dispatch on the mapping type
# (`mapping_type`, `required_geo_diff_order`, `apply_mapping!`), so the
# mapping can live here without upstream changes. The gradient rule follows
# from ∇ₓ(det J) = det J ⋅ J⁻ᵀ ⋅ (J⁻ᵀ ⊡ H) with H the geometry Hessian
# (which is why one extra geometric derivative order is required, exactly as
# for the other Piola mappings):
#
#     ∇ₓN = (∇̂N̂ ⋅ J⁻¹) / det J - N̂ ⋅ (J⁻ᵀ ⊡ (H ⋅ J⁻¹)) / det J.
#
# Only rdim == sdim cells are supported (det J), like the other Piola
# mappings. Second-order gradients (FEValues DiffOrder 2) are not
# implemented.

using Tensors: Tensors, ⊡, ⋅

"""
    L2PiolaMapping()

The L2 Piola ("integral moment preserving") mapping for scalar
interpolations: reference shape values are pushed forward as
`N(x) = N̂(ξ) / det J`. Used by declaring
`Ferrite.mapping_type(::MyInterpolation) = L2PiolaMapping()`.
"""
struct L2PiolaMapping end

Ferrite.required_geo_diff_order(::L2PiolaMapping, fun_diff_order::Int) = 1 + fun_diff_order

@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{0}, ::L2PiolaMapping, q_point::Int, mapping_values, cell
    )
    detJ = Tensors.det(Ferrite.getjacobian(mapping_values))
    @inbounds for j in 1:Ferrite.getnbasefunctions(funvals)
        funvals.Nx[j, q_point] = funvals.Nξ[j, q_point] / detJ
    end
    return nothing
end

@inline function Ferrite.apply_mapping!(
        funvals::Ferrite.FunctionValues{1}, ::L2PiolaMapping, q_point::Int, mapping_values, cell
    )
    J = Ferrite.getjacobian(mapping_values)
    H = Ferrite.gethessian(mapping_values)
    Jinv = inv(J)
    detJ = Tensors.det(J)
    # A2 = ∇ₓ(det J) / det J², cf. Ferrite's ContravariantPiolaMapping
    A2 = (Jinv' ⊡ (H ⋅ Jinv)) / detJ
    @inbounds for j in 1:Ferrite.getnbasefunctions(funvals)
        Nξ = funvals.Nξ[j, q_point]
        funvals.Nx[j, q_point] = Nξ / detJ
        funvals.dNdx[j, q_point] = (funvals.dNdξ[j, q_point] ⋅ Jinv) / detJ - Nξ * A2
    end
    return nothing
end
