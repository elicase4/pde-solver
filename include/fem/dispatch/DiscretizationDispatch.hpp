#ifndef PDESOLVER_FEM_DISPATCH_DISCRETIZATIONDISPATCH_HPP
#define PDESOLVER_FEM_DISPATCH_DISCRETIZATIONDISPATCH_HPP

#include <initializer_list>
#include <utility>

#include "core/Types.hpp"

#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"

#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadratureHex.hpp"

namespace pdesolver {
	namespace fem {
		namespace dispatch {

			enum class ElementFamily {
				Quad,
				Tri,
				Hex,
				Tet,
				Wedge
			};

			// TODO: add element tag in mesh to match the element family types
			inline ElementFamily inferElementFamily(Index npd) {
				return (npd == 2) ? ElementFamily::Quad : ElementFamily::Hex;
			}

			template<Index... Candidates, typename Cont>
			bool selectValue(Index value, Cont&& cont) {

				bool matched = false;

				(void) std::initializer_list<int>{(Candidates == value ? (cont.template operator()<Candidates>(), matched = true, 0) : 0)...};

				return matched;

			}

			template<typename Cont>
			bool selectOrder(Index p, Cont&& cont) {
				return selectValue<1, 2, 3>(p, std::forward<Cont>(cont));
			}

			template<typename Cont>
			bool selectQuadraturePoints(Index n, Cont&& cont) {
				return selectValue<1, 2, 3, 4>(n, std::forward<Cont>(cont));
			}

			template<Index NSD, typename Visitor>
			bool dispatchQuad(Index px, Index py, Index xi, Index eta, Visitor&& visitor) {

				return selectOrder(px, [&]<Index Px>() {
					return selectOrder(py, [&]<Index Py>() {
						return selectQuadraturePoints(xi, [&]<Index Xi>() {
							return selectQuadraturePoints(eta, [&]<Index Eta>() {

								using BasisT = basis::LagrangeQuad<Px, Py>;
								using QuadVolT = quadrature::GaussQuadratureQuad<Xi, Eta>;

								constexpr Index NBdy = (Xi > Eta) ? Xi : Eta;
								using QuadBdyT = quadrature::GaussQuadrature1D<NBdy>;

								return visitor.template operator()<NSD, BasisT, QuadVolT, QuadBdyT>();

							});
						});
					});
				});

			}

			template<Index NSD, typename Visitor>
			bool dispatchHex(Index px, Index py, Index pz, Index xi, Index eta, Index zeta, Visitor&& visitor) {

				return selectOrder(px, [&]<Index Px>() {
					return selectOrder(py, [&]<Index Py>() {
						return selectOrder(pz, [&]<Index Pz>() {
							return selectQuadraturePoints(xi, [&]<Index Xi>() {
								return selectQuadraturePoints(eta, [&]<Index Eta>() {
									return selectQuadraturePoints(zeta, [&]<Index Zeta>() {

										using BasisT = basis::LagrangeHex<Px, Py, Pz>;
										using QuadVolT = quadrature::GaussQuadratureHex<Xi, Eta, Zeta>;

										constexpr Index NBdy = (Xi > Eta) ? ((Xi > Zeta) ? Xi : Zeta) : ((Eta > Zeta) ? Eta : Zeta);
										using QuadBdyT = quadrature::GaussQuadratureQuad<NBdy, NBdy>;

										return visitor.template operator()<NSD, BasisT, QuadVolT, QuadBdyT>();

									});
								});
							});
						});
					});
				});

			}

			template<Index NSD, typename Visitor>
			bool dispatchNPD2(ElementFamily family, Index px, Index py, Index xi, Index eta, Visitor&& visitor) {

				if (family == ElementFamily::Quad) {
					return dispatchQuad<NSD>(px, py, xi, eta, std::forward<Visitor>(visitor));
				}

				return false;

			}

			template<Index NSD, typename Visitor>
			bool dispatchNPD3(ElementFamily family, Index px, Index py, Index pz, Index xi, Index eta, Index zeta, Visitor&& visitor) {

				if (family == ElementFamily::Hex) {
					return dispatchHex<NSD>(px, py, pz, xi, eta, zeta, std::forward<Visitor>(visitor));
				}

				return false;

			}

			template<typename Visitor>
			bool dispatch(Index nsd, Index npd, ElementFamily family, Index px, Index py, Index pz, Index xi, Index eta, Index zeta, Visitor&& visitor) {

				if (npd == 2) {

					if (nsd == 2) {
						return dispatchNPD2<2>(family, px, py, xi, eta, std::forward<Visitor>(visitor));
					}
					if (nsd == 3) {
						return dispatchNPD2<3>(family, px, py, xi, eta, std::forward<Visitor>(visitor));
					}
					return false;

				}

				if (npd == 3) {

					if (nsd == 3) {
						return dispatchNPD3<3>(family, px, py, pz, xi, eta, zeta, std::forward<Visitor>(visitor));
					}
					return false;

				}

				return false;

			}

		} // namespace dispatch
	} // namespace fem
} // namespace pdesolver

#endif
