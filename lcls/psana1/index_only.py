from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from dials.algorithms.indexing import index_reflections
from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
import sys

expts = ExperimentListFactory.from_json_file(sys.argv[1], check_format=False)
refls = flex.reflection_table.from_file(sys.argv[2])
outname = sys.argv[3]

indexed = flex.reflection_table()
for expt_id, experiment in enumerate(expts):
    subset = refls.select(refls['id'] == expt_id)
    subset["id"] = flex.int(len(subset), -1)
    subset["imageset_id"] = flex.int(len(subset), 0)
    subset.centroid_px_to_mm(expts[expt_id:expt_id+1])
    subset.map_centroids_to_reciprocal_space(expts[expt_id:expt_id+1])
    index_reflections(subset, expts[expt_id:expt_id+1])
    subset = subset.select(subset['id'] >= 0)
    ref_predictor = ExperimentsPredictorFactory.from_experiments(expts[expt_id:expt_id+1])
    ref_predictor(subset)
    subset['id'] = flex.int(len(subset), expt_id)
    indexed.extend(subset)
indexed["entering"] = flex.bool(len(indexed), False)
indexed.as_msgpack_file(outname)
