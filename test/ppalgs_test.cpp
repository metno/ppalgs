#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ppalgs test
#include <boost/test/unit_test.hpp>

#define BOOST_TEST_LOG_LEVEL "all"

#include "Ducting.h"
#include "Icing.h"
#include "Contrails.h"
#include "GribHandler.h"
#include "NetCDFHandler.h"

struct DuctingFixture {
	DuctingFixture() { BOOST_TEST_MESSAGE( "setup fixture" ); }
    ~DuctingFixture() { BOOST_TEST_MESSAGE( "teardown fixture" ); }

    ///XXX
    //Ducting ducting;
};

struct IcingFixture {
	IcingFixture() { BOOST_TEST_MESSAGE( "setup fixture" ); }
    ~IcingFixture() { BOOST_TEST_MESSAGE( "teardown fixture" ); }

    ///XXX
    //Icing icing;
};

struct ContrailsFixture {
	ContrailsFixture() { BOOST_TEST_MESSAGE( "setup fixture" ); }
    ~ContrailsFixture() { BOOST_TEST_MESSAGE( "teardown fixture" ); }

    ///XXX
    //Contrails contrails;
};

BOOST_FIXTURE_TEST_SUITE(TestDucting, DuctingFixture)

BOOST_AUTO_TEST_CASE(TestDuctingStub)
{
	// no tests implemented yet

	/* use BOOST_CHECK_EQUAL_COLLECTIONS?
	int col1 [] = { 1, 2, 3, 4, 5, 6, 7 };
	int col2 [] = { 1, 2, 4, 4, 5, 7, 7 };
	BOOST_CHECK_EQUAL_COLLECTIONS( col1, col1+7, col2, col2+7 );
	*/
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(TestIcing, IcingFixture)

BOOST_AUTO_TEST_CASE(TestIcingStub)
{
	// no tests implemented yet

	/* use BOOST_CHECK_EQUAL_COLLECTIONS?
	int col1 [] = { 1, 2, 3, 4, 5, 6, 7 };
	int col2 [] = { 1, 2, 4, 4, 5, 7, 7 };
	BOOST_CHECK_EQUAL_COLLECTIONS( col1, col1+7, col2, col2+7 );
	*/
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(TestContrails, ContrailsFixture)

BOOST_AUTO_TEST_CASE(TestContrailsStub)
{
	// no tests implemented yet

	/* use BOOST_CHECK_EQUAL_COLLECTIONS?
	int col1 [] = { 1, 2, 3, 4, 5, 6, 7 };
	int col2 [] = { 1, 2, 4, 4, 5, 7, 7 };
	BOOST_CHECK_EQUAL_COLLECTIONS( col1, col1+7, col2, col2+7 );
	*/
}

BOOST_AUTO_TEST_SUITE_END()
