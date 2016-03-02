%define name LLL
%define version 0.9
%define unmangled_version 0.9
%define release 1

Summary: LLL - libraries for Laue lens calculation
Name: %{name}
Version: %{version}
Release: %{release}
Source0: %{name}-%{unmangled_version}.tar.gz
License: UNKNOWN
Group: Development/Libraries
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
Prefix: %{_prefix}
BuildArch: noarch
Vendor: Alessandro Pisa <alessandro...pisa@@@gmail...com>
Url: http://darkmoon.altervista.org

%description

LLL - libraries for Laue lens calculation


%prep
%setup -n %{name}-%{unmangled_version}

%build
python setup.py build

%install
python setup.py install --root=$RPM_BUILD_ROOT --record=INSTALLED_FILES

%clean
rm -rf $RPM_BUILD_ROOT

%files -f INSTALLED_FILES
%defattr(-,root,root)
