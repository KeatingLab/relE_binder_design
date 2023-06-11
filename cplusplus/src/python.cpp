#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "mstfuser.h"
#include "msttypes.h"
#include "mstoptions.h"

#include "bridgeseeds.h"
#include "generateseeds.h"
#include "residuecontact.h"
#include "utilities.h"

// Source: https://stackoverflow.com/questions/15842126/feeding-a-python-list-into-a-function-taking-in-a-vector-with-boost-python/15940413#15940413

/// @brief Type that allows for registration of conversions from
///                python iterable types.
struct iterable_converter
{
    /// @note Registers converter from a python interable type to the
    ///             provided type.
    template <typename Container>
    iterable_converter&
    from_python() {
        boost::python::converter::registry::push_back(
                                                      &iterable_converter::convertible,
                                                      &iterable_converter::construct<Container>,
                                                      boost::python::type_id<Container>());

        // Support chaining.
        return *this;
    }

    /// @brief Check if PyObject is iterable.
    static void* convertible(PyObject* object) {
        return PyObject_GetIter(object) ? object : NULL;
    }

    /// @brief Convert iterable PyObject to C++ container type.
    ///
    /// Container Concept requirements:
    ///
    ///     * Container::value_type is CopyConstructable.
    ///     * Container can be constructed and populated with two iterators.
    ///         I.e. Container(begin, end)
    template <typename Container>
    static void construct(
                          PyObject* object,
                          boost::python::converter::rvalue_from_python_stage1_data* data) {
        namespace python = boost::python;
        // Object is a borrowed reference, so create a handle indicting it is
        // borrowed for proper reference counting.
        python::handle<> handle(python::borrowed(object));

        // Obtain a handle to the memory block that the converter has allocated
        // for the C++ type.
        typedef python::converter::rvalue_from_python_storage<Container>
        storage_type;
        void* storage = reinterpret_cast<storage_type*>(data)->storage.bytes;

        typedef python::stl_input_iterator<typename Container::value_type>
        iterator;

        // Allocate the C++ type into the converter's memory block, and assign
        // its handle to the converter's convertible variable.    The C++
        // container is populated by passing the begin and end iterators of
        // the python object to the container's constructor.
        new (storage) Container(
                                iterator(python::object(handle)), // begin
                                iterator()); // end
        data->convertible = storage;
    }
};

vector<Residue*> getSeedResWithResection(Structure* seed, bool nterminus = true, int offset = 0) {
    vector<Residue*> all_res = seed->getResidues();
    vector<Residue*> ret_res;
    for (int i = 0; i < all_res.size(); i++) {
        if (nterminus && (i < offset)) {
            // skip n-terminal residue
        } else if (!nterminus && (i >= (seed->residueSize() - offset))) {
            // skip c-terminal residue
        } else {
            ret_res.push_back(all_res[i]);
        }
    }
    return ret_res;
}

/*
 * This is a macro Boost.Python provides to signify a Python extension module.
 */
BOOST_PYTHON_MODULE(peptide_binder_design) {
    // An established convention for using boost.python.
    using namespace boost::python;

   iterable_converter()
    .from_python<vector<Atom*>>()
    .from_python<vector<Residue*>>()
    .from_python<vector<Structure*>>()
    .from_python<vector<shared_ptr<Structure>>>()
    // .from_python<vector<double>>()
    // .from_python<vector<int>>()
    // .from_python<vector<vector<int>>>()
    // .from_python<vector<vector<Residue*>>>()
    // .from_python<vector<std::string>>()
    // .from_python<vector<MST::Sequence>>()
    ;
    // Various translations of vectors to lists needed by the classes below
    class_<vector<Atom*>>("AtomList")
        .def(vector_indexing_suite<vector<Atom*>>());

    class_<vector<Residue*>>("ResidueList")
        .def(vector_indexing_suite<vector<Residue*>>());

    class_<vector<Structure*>>("StructureList")
        .def(vector_indexing_suite<vector<Structure*>>());

    class_<vector<shared_ptr<Structure>>>("spStructureList")
        .def(vector_indexing_suite<shared_ptr<Structure>>());

    // class_<vector<int>>("IntList")
    //     .def(vector_indexing_suite<vector<int>>());

    // class_<vector<mstreal>>("PointList")
    //     .def(vector_indexing_suite<vector<double>>());

    // class_<vector<std::string>>("StringList")
    //     .def(vector_indexing_suite<vector<std::string>>());

    // class_<vector<Sequence>>("SequenceList")
    //     .def(vector_indexing_suite<vector<MST::Sequence>>());

    class_<MST::Residue>("Residue", init<>())
        .def("atomSize", &MST::Residue::atomSize)
        .add_property("name", &MST::Residue::getName)
        .add_property("num", &MST::Residue::getNum)
        .def("getAtom", &MST::Residue::getAtom, return_value_policy<reference_existing_object>())
        .def("previousResidue", &MST::Residue::previousResidue, return_value_policy<reference_existing_object>())
        .def("nextResidue", &MST::Residue::nextResidue, return_value_policy<reference_existing_object>())
        .add_property("phi", &MST::Residue::getPhi)
        .add_property("psi", &MST::Residue::getPsi)
        .add_property("omega", &MST::Residue::getOmega)
        .def("isBadDihedral", &MST::Residue::isBadDihedral)
        .staticmethod("isBadDihedral")
        .def("areBonded", static_cast<bool (*) (const MST::Residue&, const MST::Residue&, MST::mstreal)>(&MST::Residue::areBonded))
        .staticmethod("areBonded")
        .def("getResidueIndex", &MST::Residue::getResidueIndex)
        .def("getResidueIndexInChain", &MST::Residue::getResidueIndexInChain)
        .def("getStructure", &MST::Residue::getStructure, return_value_policy<reference_existing_object>())
        .def("getParent", &MST::Residue::getParent, return_value_policy<reference_existing_object>())
    ;

    class_<Structure>("Structure",
                      init<string,string>())
        .def(init<vector<Residue*>>())
        .def("chainSize", &MST::Structure::chainSize)
        .def("residueSize", &MST::Structure::residueSize)
        .def("atomSize", &MST::Structure::atomSize)
        .def("getChain", &MST::Structure::getChain, return_value_policy<reference_existing_object>())
        .def("getResidue", &MST::Structure::getResidue,
            return_value_policy<reference_existing_object>())
        .def("getAtoms", &MST::Structure::getAtoms)
        .def("getResidues", &MST::Structure::getResidues)
        .def("appendChain", +[](MST::Structure& structure, MST::Chain *chain) {
            structure.appendChain(new MST::Chain(*chain));
        })/*static_cast<bool (MST::Structure::*) (Chain*, bool)>(&MST::Structure::appendChain))*/
        .def("deleteChain", &MST::Structure::deleteChain)
        .def("addAtom", static_cast<void (MST::Structure::*) (MST::Atom*)>(&MST::Structure::addAtom))
        .def("addAtoms", static_cast<void (MST::Structure::*) (std::vector<MST::Atom*>*)>(&MST::Structure::addAtoms))
        .def("addResidue", &MST::Structure::addResidue, return_value_policy<manage_new_object>())
        .def("getResidueIndex", &MST::Structure::getResidueIndex)
        .def("__eq__", &MST::Structure::operator==)
        .def("__ne__", &MST::Structure::operator!=)
        // the static_cast is needed to disambiguate an overloaded function
        .def("writePDB", static_cast<void (MST::Structure::*) (const std::string&, std::string) const>(&MST::Structure::writePDB))
        // .def("__str__", &Py_Structure::structureToString)
        .add_property("name", &MST::Structure::getName, &MST::Structure::setName)
        .def("reassignChainsByConnectivity", static_cast<MST::Structure (MST::Structure::*) (MST::mstreal)> (&MST::Structure::reassignChainsByConnectivity))
    ;

    class_<seedBinaryFile>("seedBinaryFile",init<string,bool,bool>())
        .def("next",&seedBinaryFile::next,
            return_internal_reference<>())
        .def("hasNext",&seedBinaryFile::hasNext)
        .def("skip",&seedBinaryFile::skip)
        .def("reset",&seedBinaryFile::reset)
        .def("getStructureNamed",&seedBinaryFile::getStructureNamed
            , return_internal_reference<>())
    ;

    class_<clashChecker>("clashChecker")
        .def("setStructure",&clashChecker::setStructure)
        .def("setResidues",&clashChecker::setResidues)
        .def("checkForClashesToStructure",&clashChecker::checkForClashesToStructure)
    ;

    class_<findSeedBridge>("findSeedBridge",init<string,string,int>())
        .def("setMaxSeedDistance",&findSeedBridge::setMaxSeedDistance)
        .def("setSearchQuery",&findSeedBridge::setSearchQuery)
        .def("searchSeedsByCADistance",&findSeedBridge::searchSeedsByCADistance)
        .def("loadStructureDB",&findSeedBridge::loadStructureDB)
        .def("verifyMatchesBySuperPosition",&findSeedBridge::verifyMatchesBySuperposition)
        .def("reportAPVBoundaries",&findSeedBridge::reportAPVBoundaries)
        .def("getVerifiedBridgeStructures",&findSeedBridge::getVerifiedBridgeStructures,return_value_policy<manage_new_object>())
    ;

    class_<fuseSeedsAndBridge>("fuseSeedsAndBridge",init<int,string>())
        .def("setClashCheck",&fuseSeedsAndBridge::setClashCheck)
        .def("setSeeds",&fuseSeedsAndBridge::setSeeds)
        .def("setBridgeStructures",&fuseSeedsAndBridge::setBridgeStructures)
        .def("fuse",&fuseSeedsAndBridge::fuse,return_value_policy<manage_new_object>())
    ;

    def("getSeedResWithResection", getSeedResWithResection);
}